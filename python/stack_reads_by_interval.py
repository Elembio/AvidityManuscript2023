import os
import gc
import sys
import gzip
import pysam
import random
import logging
import argparse
import numpy as np
from tqdm import tqdm

from util import get_called_and_aligned, process_sequence_records, NUCLEOTIDE2IDX


MAPQ_THRESH = 0

logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)
log = logging.getLogger("stack-reads-by-interval")


def get_overlap(a, b):
    """test to for overlap between two intervals.
    """

    if not a[0] < a[1]:
        log.warning(f"invalid interval a - skipping ({a[0]} > {a[1]})")
        return None

    if not b[0] < b[1]:
        log.warning(f"invalid interval b - skipping ({b[0]} > {b[1]})")
        return None
    
    assert a[0] < a[1], f"{a[0]} < {a[1]}"
    assert b[0] < b[1], f"{b[0]} < {b[1]}"

    if a[0] <= b[0] and a[1] > b[0]:
        return b[0], min(a[1], b[1])
    
    if b[0] <= a[0] and b[1] > a[0]:   
        return a[0], min(a[1], b[1])
           
    return None

def trim_read_idx(read_iter, num_reads):
    """
    trim readidx2intervalidx arr to match mats
    """

    result = np.zeros((num_reads), dtype=np.uint32)
    
    for read_idx, record in enumerate(read_iter):
        if record is None:
            continue
        result[read_idx] = record

    return result


def main():
    
    parser = argparse.ArgumentParser(description='stack-reads-by-interval')

    parser.add_argument('--bam', type=str, required=True, default=None, help='bam file - pos sorted + subset to only overlapping read')
    parser.add_argument('--bed', type=str, required=True, default=None, help='bed file - sorted by *.genome')
    parser.add_argument('-o', '--output-prefix', type=str, required=True, default=None, help='output prefix - workflow name')

    parser.add_argument('--chr', type=str, default=None, help='optional chr by which to chunk input')
    parser.add_argument('--rate', type=float, default=1.0, help='rate at which to select/retain reads')
    parser.add_argument('--slop', type=int, default=5, help='BED interval slop [5]')
    parser.add_argument('--min-interval', type=int, default=1, help='Minimum num bases observed of bed interval [1]')
    parser.add_argument('--min-pre', type=int, default=10, help='Minimum num bases observed after bed interval [10]')
    parser.add_argument('--min-post', type=int, default=10, help='Minimum num bases observed before bed interval [10]')
    parser.add_argument('--pre-context', type=int, default=None, help='Size of context to track in pre-interval sequence [Inf]')
    parser.add_argument('--post-context', type=int, default=None, help='Size of context to track in post-interval sequence [Inf]')
    parser.add_argument('--window', type=int, default=150, help='Window side for pre and post homopolymer context [150]')
    parser.add_argument('--read_length', type=int, default=150, help='Read Length')
    parser.add_argument('--homopolymer', action='store_true', help='enable homopolymer mode')
    parser.add_argument('--debug', action='store_true', help='enable debug mode')

    args = parser.parse_args()

    log.setLevel(logging.INFO)
    if args.debug:
        log.setLevel(logging.DEBUG)

    args.window = min(args.window,args.window-args.min_interval-max(args.min_pre,args.min_post))

    print(args.pre_context)

    if not args.pre_context:
        args.pre_context = args.read_length - args.min_interval - args.min_pre
    if not args.post_context:
        args.post_context = args.read_length - args.min_interval - args.min_post

    log.info(f"bam: {args.bam}")
    log.info(f"bed: {args.bed}")
    log.info(f"chr: {args.chr}")
    log.info(f"output_prefix: {args.output_prefix}")
    log.info(f"rate: {args.rate}")
    log.info(f"slop: {args.slop}")
    log.info(f"min-interval: {args.min_interval}")
    log.info(f"min-pre: {args.min_pre}")
    log.info(f"min-post: {args.min_post}")
    log.info(f"pre-context: {args.pre_context}")
    log.info(f"post-context: {args.post_context}")
    log.info(f"window: {args.window}")
    log.info(f"read_length: {args.read_length}")
    log.info(f"homopolymer: {args.homopolymer}")

    contig2index=dict()

    with pysam.Samfile(args.bam, "rb") as bam_fh:
        # compute total reads
        total_reads = 0
        for i,idxstat in enumerate(bam_fh.get_index_statistics()):
            contig = idxstat.contig
            contig2index[contig]=i
            if args.chr and contig != args.chr:
                continue
            
            # only count reads for chr of interest
            total_reads += idxstat.total
    
    total_reads = total_reads
    chunk_size = max(1,int(total_reads/100))
    log.info(f"total_reads: {total_reads}")
    log.info(f"chunk_size: {chunk_size}")

    # open bed file

    # Preallocate arrays
    pre_called_sequences = [None] * total_reads
    pre_aligned_sequences = [None] * total_reads
    post_called_sequences = [None] * total_reads
    post_aligned_sequences = [None] * total_reads
    homoP_bases = [None] * total_reads
    interval_lengths = [None] * total_reads
    read_to_interval_idx = [None] * total_reads
    read_indexes = [None] * total_reads
    read_orientations = [None] * total_reads
    
    # Index into above arrays
    read_idx = 0
    interval_idx = -1

    # Track max length of pre and post homo-P sequences
    # Deletions are filled with 'D' so can be longer than original read length
    max_pre_seq_len = 0
    max_post_seq_len = 0

    bed_intervals = []

    with gzip.open(args.bed, 'rt') as in_bed_fh:
        
        #get first bed line
        bed_line = in_bed_fh.readline().rstrip()

        # skip over comment lines
        while(bed_line.startswith("#")):
            bed_line = in_bed_fh.readline()
            if not bed_line:
                break
        
        # get initial bed line
        bed_chr,bed_start,bed_end = bed_line.rstrip().split("\t")
        log.info(f"initial {interval_idx} bed - {bed_chr}:{bed_start}-{bed_end}")
        bed_start, bed_end = int(bed_start) + args.slop, int(bed_end) - args.slop
        log.info(f"initial {interval_idx} bed w/o slop - {bed_chr}:{bed_start}-{bed_end}")
        bed_intervals.append((bed_chr, bed_start, bed_end))
        interval_idx += 1 

        # open bam file
        bam_fh = pysam.Samfile(args.bam, "rb")    
        num_secondary, num_supplementary, num_unmapped, num_mapped, num_below_mapq_thresh, num_rate_skipped = 0,0,0,0,0,0
        num_skipped_interval_len, num_skipped_pre_interval_len, num_skipped_post_interval_len, num_slop_overlap_misses = 0,0,0,0
    
        # read bam 
        read_overlaps = 0
        log.info("reading bam ... ")
        with tqdm(total=total_reads) as pbar:
            
            for i,read in enumerate(bam_fh.fetch(args.chr)):
                
                # update progress bar
                if i % chunk_size == 0:
                    if not args.debug:
                        pbar.update(chunk_size)
                
                skip=False

                # skip non-unique mapped reads
                if read.is_secondary == True:
                    num_secondary+=1
                    skip=True
                if read.is_supplementary == True:
                    num_supplementary+=1
                    skip=True
                if read.is_unmapped == True:
                    num_unmapped+=1
                    skip=True
                if read.mapping_quality <= MAPQ_THRESH:
                    num_below_mapq_thresh+=1
                    skip=True
                
                if skip:
                    continue

                num_mapped += 1

                query_name = read.query_name
                bam_chr = bam_fh.get_reference_name(read.reference_id)
                bam_start = read.reference_start
            
                # find overlapping bed interval
                while contig2index[bed_chr] < contig2index[bam_chr]:
                    bed_line = in_bed_fh.readline()
                    if not bed_line:
                        break
                    bed_chr,bed_start,bed_end = bed_line.rstrip().split("\t")
                    bed_start, bed_end = int(bed_start) + args.slop, int(bed_end) - args.slop
                    bed_intervals.append(None)
                    interval_idx += 1
                
                #
                # we should now be on the same chr between bam/bed
                #

                called_seq, aligned_seq = get_called_and_aligned(read, skip_ins=False)
                
                # in case of deletions, get read range here
                # deletions are filled in with 'D' in the called sequence
                bam_end = bam_start + len(called_seq)

                # Read is skipped if cigar string contains operations other than CMATCH, CSOFT_CLIP, CINS, CDEL
                if not called_seq:
                    log.debug("Cigar string doesn't meet requirements, skipping read")
                    continue
                    
                # handle bam edge(s)
                if contig2index[bam_chr] < contig2index[bed_chr]:
                    log.debug(f"skipping bam-edge - i:{i} chunk:{int(i/chunk_size)} read_idx:{read_idx} - {query_name}\t{bam_chr}:{bam_start}-{bam_end}\t{bed_chr}:{bed_start}-{bed_end}")
                    continue

                assert bed_chr == bam_chr
                
                while bed_chr == bam_chr and bed_end < bam_start:
                    bed_line = in_bed_fh.readline()
                    if not bed_line:
                        break
                    bed_chr,bed_start,bed_end = bed_line.rstrip().split("\t")
                    bed_start, bed_end = int(bed_start) + args.slop, int(bed_end) - args.slop
                    bed_intervals.append((bed_chr, bed_start, bed_end))
                    interval_idx += 1
                    
                if bed_chr != bam_chr:
                    continue
                
                log.debug(f"line:{i} - chunk-{int(i/chunk_size)} - reads:[{read_idx}] - interval:[{interval_idx}] bam={bam_chr}:{bam_start}-{bam_end} x bed={bed_chr}:{bed_start}-{bed_end}")

                overlap = get_overlap([bam_start,bam_end], [int(bed_start),int(bed_end)])

                if not overlap:
                    num_slop_overlap_misses += 1
                    continue
                
                # only continue if random < rate
                rand = random.random()
                if rand > args.rate:
                    num_rate_skipped+=1 
                    continue

                # Get lengths of three contigs:
                # pre_homop := sequence preceding homopolymer
                # homop := homopolymer sequence indicated by bed interval
                # post_homop := sequence succeeding homopolymer
                if read.is_reverse:
                    pre_interval_len = max(0, bam_end - overlap[1])
                    post_interval_len = max(0, overlap[0] - bam_start)
                else:
                    pre_interval_len = max(0, overlap[0] - bam_start)
                    post_interval_len =  max(0, bam_end - overlap[1])
                
                interval_len = len(called_seq) - (pre_interval_len + post_interval_len)

                if args.min_interval > interval_len or args.min_pre > pre_interval_len or args.min_post > post_interval_len:
                    if interval_len < args.min_interval:
                        num_skipped_interval_len +=1
                        continue
                    if pre_interval_len < args.min_pre:
                        num_skipped_pre_interval_len +=1
                        continue
                    if post_interval_len < args.min_post:
                        num_skipped_post_interval_len +=1
                        continue

                # Update max lengths
                if pre_interval_len > max_pre_seq_len:
                    max_pre_seq_len = pre_interval_len
                if post_interval_len > max_post_seq_len:
                    max_post_seq_len = post_interval_len
                
                homoP_base = "N"
                if args.homopolymer:
                    # Find homopolymer base
                    # The base recorded is according to instrument-cycle order
                    homoP_base_set = set(aligned_seq[pre_interval_len:pre_interval_len+interval_len])
                    if 'N' in homoP_base_set:
                        homoP_base_set.remove('N')

                    # If multiple homopolymers are found in the same interval, skip this read
                    # This happens as a result of BED dedup on intervals with nonzero slop (GIAB issue)
                    if len(homoP_base_set) != 1:
                        continue
                    
                    homoP_base = homoP_base_set.pop()

                # Count this read as overlapping a valid homo-P interval
                read_overlaps += 1

                log.debug(f"line:{i} - chunk-{int(i/chunk_size)} - reads:[{read_idx}] - interval:[{interval_idx}] interval_len:{interval_len} pre_interval_len:{pre_interval_len} post_interval_len:{post_interval_len} bam={bam_chr}:{bam_start}-{bam_end} x bed={bed_chr}:{bed_start}-{bed_end}")

                # Define array start/end for pre and post context sequences
                pre_start = max(0, pre_interval_len - args.post_context)
                pre_end = pre_interval_len
                post_start = pre_interval_len + interval_len
                post_end = pre_start + args.pre_context

                # Record pre_homop and post_homop called and aligned sequences
                # The bases recorded are according to instrument-cycle order
                pre_called_sequences[read_idx] = called_seq[pre_start:pre_end]
                pre_aligned_sequences[read_idx] = aligned_seq[pre_start:pre_end]
                post_called_sequences[read_idx] = called_seq[post_start:post_end]
                post_aligned_sequences[read_idx] = aligned_seq[post_start:post_end]

                # Metadata
                homoP_bases[read_idx] = homoP_base
                interval_lengths[read_idx] = interval_len
                read_to_interval_idx[read_idx] = interval_idx
                read_indexes[read_idx] = i
                read_orientations[read_idx] = read.is_reverse

                read_idx += 1

        bam_fh.close()
                
        log.info(f"total_reads: {total_reads}")
        log.info(f"num_secondary: {num_secondary}")
        log.info(f"num_supplementary: {num_supplementary}")
        log.info(f"num_unmapped: {num_unmapped}")
        log.info(f"num_below_mapq_{MAPQ_THRESH}: {num_below_mapq_thresh}")
        log.info(f"num_mapped: {num_mapped}")
        log.info(f"num_rate_skipped: {num_rate_skipped}")
        log.info(f"num_slop_overlap_misses: {num_slop_overlap_misses}")
        log.info(f"num_skipped_interval_len: {num_skipped_interval_len}")
        log.info(f"num_skipped_pre_interval_len: {num_skipped_pre_interval_len}")
        log.info(f"num_skipped_post_interval_len: {num_skipped_post_interval_len}")
        
        log.info(f"read_overlaps: {read_overlaps}")
        log.info(f"read_idx: {read_idx}")
        log.info(f"pre-called stack height: {sum(x is not None for x in pre_called_sequences)}")
        log.info(f"pre-aligned stack height: {sum(x is not None for x in pre_aligned_sequences)}")
        log.info(f"post-called stack height: {sum(x is not None for x in post_called_sequences)}")
        log.info(f"post-aligned stack height: {sum(x is not None for x in post_aligned_sequences)}")

        log.info("Building base/mismatch tsv...")
        
        # right_shift = True creates offset so that position -1 is the last basecall in pre_homop sequence
        pre_called_mat = process_sequence_records(pre_called_sequences, max_pre_seq_len, read_idx, right_shift = True)
        del pre_called_sequences
        pre_aligned_mat = process_sequence_records(pre_aligned_sequences, max_pre_seq_len, read_idx, right_shift = True)
        del pre_aligned_sequences
        gc.collect()

        pre_mismatch_mat, pre_base_mask = get_mismatch_matrix(pre_called_mat, pre_aligned_mat)
        del pre_called_mat
        del pre_aligned_mat
        gc.collect()

        post_called_mat = process_sequence_records(post_called_sequences, max_post_seq_len, read_idx)
        del post_called_sequences
        post_aligned_mat = process_sequence_records(post_aligned_sequences, max_post_seq_len, read_idx)
        del post_aligned_sequences
        gc.collect()

        post_mismatch_mat, post_base_mask = get_mismatch_matrix(post_called_mat, post_aligned_mat)
        del post_called_mat
        del post_aligned_mat
        gc.collect()

        usable_read_to_interval_idx = trim_read_idx(read_to_interval_idx, read_idx)
        
        if len(pre_base_mask) == 0 or len(post_base_mask) == 0:
            log.error(f"no data found for chromosome {args.chr}")
            sys.exit()

        pre_error_ct = pre_mismatch_mat.sum()
        pre_correct_ct = len(pre_mismatch_mat) - pre_error_ct
        pre_base_ct = pre_base_mask.sum()
        pre_error_rate = (pre_error_ct / pre_base_ct) * 100
        
        post_error_ct = post_mismatch_mat.sum()
        post_correct_ct = len(post_mismatch_mat) - post_error_ct
        post_base_ct = post_base_mask.sum()
        post_error_rate = (post_error_ct / post_base_ct) * 100
        
        pre_error_by_pos = np.zeros(max_pre_seq_len)
        pre_weight_by_pos = np.zeros_like(pre_error_by_pos)

        post_error_by_pos = np.zeros(max_post_seq_len)
        post_weight_by_pos = np.zeros_like(post_error_by_pos)

        # dump error by offset
        log.info("dumping offset-error")
        offset_error_tsv = os.path.join(args.output_prefix+".offset-error.tsv")
        with open(offset_error_tsv, "w") as out_tsv:

            for pos in range(max_pre_seq_len):
                pre_idx = pos - max_pre_seq_len
                pre_weight_by_pos[pre_idx] = pre_base_mask[:,pre_idx].sum()
                pre_error_by_pos[pre_idx] = (pre_mismatch_mat[:,pre_idx].sum() / pre_weight_by_pos[pre_idx]) * 100
                print(args.chr,pre_idx,pre_mismatch_mat[:,pre_idx].sum(),pre_weight_by_pos[pre_idx],pre_error_by_pos[pre_idx],sep="\t",file=out_tsv)

            for pos in range(max_post_seq_len):                
                pos_idx = pos
                post_weight_by_pos[pos_idx] = post_base_mask[:,pos_idx].sum()
                post_error_by_pos[pos_idx] = (post_mismatch_mat[:,pos_idx].sum() / post_weight_by_pos[pos_idx]) * 100
                print(args.chr,pos_idx,post_mismatch_mat[:,pos_idx].sum(),post_weight_by_pos[pos_idx],post_error_by_pos[pos_idx],sep="\t",file=out_tsv)
        
        log.info("completed offset-error")


        log.info("dumping interval-error")
        interval_error_tsv = os.path.join(args.output_prefix+".interval-error.tsv")
        with open(interval_error_tsv, "w") as out_tsv:
            unique_intervals = np.unique(usable_read_to_interval_idx)
            log.info(f"num_intervals: {len(unique_intervals)}")
            
            for i,interval in enumerate(unique_intervals[:-1]):
                interval += 1 
    
                interval_mask = usable_read_to_interval_idx == interval
                
                total_before_error = pre_mismatch_mat[interval_mask].sum() 
                total_after_error = post_mismatch_mat[interval_mask].sum()
                total_before_base = pre_base_mask[interval_mask].sum() 
                total_after_base = post_base_mask[interval_mask].sum()

                error_rate_before = 0.
                if total_before_base > 0:
                    error_rate_before = (total_before_error / total_before_base) * 100
                
                error_rate_after = 0.
                if total_after_base > 0:
                    error_rate_after = (total_after_error / total_after_base) * 100
    
                before_coverage = pre_base_mask[interval_mask].shape[0]
                after_coverage = post_base_mask[interval_mask].shape[0]
                
                interval_coord = f"{bed_intervals[interval][0]}:{bed_intervals[interval][1]}-{bed_intervals[interval][2]}"
                
                print(interval,interval_coord,before_coverage,total_before_error,total_before_base,error_rate_before,after_coverage,total_after_error,total_after_base,error_rate_after,sep="\t",file=out_tsv)

        log.info("completed interval-error")
        log.info(f'{args.output_prefix}\t{args.chr}\tpre-interval\t{pre_correct_ct}\t{pre_error_ct}\t{pre_error_rate}')
        log.info(f'{args.output_prefix}\t{args.chr}\tpost-interval\t{post_correct_ct}\t{post_error_ct}\t{post_error_rate}')


def get_mismatch_matrix(called_iter, aligned_iter):
    print(f"get_mismatch_matrix: called_iter:{called_iter.shape} aligned_iter:{aligned_iter.shape}")

    assert called_iter.shape == aligned_iter.shape
    no_call_mask = np.logical_and(called_iter != NUCLEOTIDE2IDX['N'], aligned_iter != NUCLEOTIDE2IDX['N'])
    del_mask = called_iter != NUCLEOTIDE2IDX['D']
    relevant_basecalls_mask = np.logical_and(no_call_mask, del_mask)

    del no_call_mask
    del del_mask
    gc.collect()

    result = np.zeros_like(called_iter, dtype=bool)
    
    print("building mismatch matrix..")
    result[relevant_basecalls_mask] = called_iter[relevant_basecalls_mask] != aligned_iter[relevant_basecalls_mask]

    return result, relevant_basecalls_mask


if __name__ == "__main__":
    main()
