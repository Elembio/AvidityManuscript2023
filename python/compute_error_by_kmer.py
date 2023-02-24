import os
import pysam
import random
import logging
import argparse
import numpy as np 
from tqdm import tqdm
from itertools import product

from util import get_called_and_aligned, byte_encode_seq, NUCLEOTIDE2IDX
from util.genome_mask import GenomeMask, SAMPLEDATA

logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)
log = logging.getLogger("compute-error-by-kmer")


def main():
    
    parser = argparse.ArgumentParser(description='compute-error-by-kmer')

    parser.add_argument('--bam', type=str, required=True, default=None, help='bam file - pos sorted + subset to only overlapping read')
    parser.add_argument('--kmer_length', type=int, default=8, help='length of kmer to compute error for')
    parser.add_argument('-o', '--output-prefix', type=str, required=True, default=None, help='output prefix - workflow name')
    parser.add_argument('--sample_name', type=str, required=False, default=None, help='sample name, coriell ID, HG001/HG002/HG005 etc')
    
    parser.add_argument('--min_base_q', type=int, required=False, default=0, help='optionally filter by min base q')
    parser.add_argument('--min_map_q', type=int, required=False, default=0, help='optionally filter by min map q')
    parser.add_argument('--chr', type=str, required=False, default=None, help='optional chr by which to chunk input')
    parser.add_argument('--rate', type=float, default=1.0, help='rate at which to select/retain reads')
    parser.add_argument('--read_length', type=int, default=150, help='Read Length')
    parser.add_argument('--debug', action='store_true', help='enable debug mode')

    args = parser.parse_args()

    log.setLevel(logging.INFO)
    if args.debug:
        log.setLevel(logging.DEBUG)

    log.info(f"bam: {args.bam}")
    log.info(f"kmer_length: {args.read_length}")
    log.info(f"output_prefix: {args.output_prefix}")
    log.info(f"chr: {args.chr}")
    log.info(f"rate: {args.rate}")
    log.info(f"read_length: {args.read_length}")
    log.info(f"sample_name: {args.sample_name}")

    if args.sample_name not in SAMPLEDATA:
        log.error(f"sample name [{args.sample_name}] not supported")
        sys.exit()

    # Build mask for skipping variants if sample name is provided
    log.info(f"loading genome mask")
    genome_mask = GenomeMask.from_sample_name(args.sample_name)
    log.info("done")

    # compute total number of reads for timing
    log.info("computing total reads")
    with pysam.Samfile(args.bam, "rb") as bam_fh:
        # compute total reads
        total_reads = 0
        for i,idxstat in enumerate(bam_fh.get_index_statistics()):
            contig = idxstat.contig
            if args.chr and contig != args.chr:
                continue
            
            # only count reads for chr of interest
            total_reads += idxstat.total
    
    # compute tqdm chunk/logging
    total_reads = total_reads
    chunk_size = max(1,int(total_reads/10000))
    log.info(f"total_reads: {total_reads}")
    log.info(f"chunk_size: {chunk_size}")

    # compute all possible kmer of length N, pre-fill dict for total/mismatch
    kmer_dict_size = 4**args.kmer_length
    if kmer_dict_size > 2**32:
          log.error(f"too many kmers {kmer_dict_size} > {2**32} limit")
          raise SystemExit(f"too many kmers {kmer_dict_size} > {2**32} limit")

    kmer_vocabulary = [''.join(c) for c in product('ACTG', repeat=args.kmer_length)]
    kmer_dict = {}
    for kmer in kmer_vocabulary:
        kmer_dict[kmer] = {}
        kmer_dict[kmer]["all"] = np.zeros(args.kmer_length, dtype=int)            
        kmer_dict[kmer]["mismatch"] = np.zeros(args.kmer_length, dtype=int)            
        kmer_dict[kmer]["correct"] = 0
        kmer_dict[kmer]["incorrect"] = 0

    assert kmer_dict_size == len(kmer_vocabulary)
    log.info(f"kmer_dict_size={kmer_dict_size}")

    # open bam file
    bam_fh = pysam.Samfile(args.bam, "rb")    
    num_secondary, num_supplementary, num_unmapped, num_mapped, num_below_mapq_thresh, num_rate_skipped, num_kmers_used = 0,0,0,0,0,0,0

    # read bam 
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
            if read.mapping_quality <= args.min_map_q:
                num_below_mapq_thresh+=1
                skip=True
            
            if skip:
                continue

            # only continue if random < rate
            rand = random.random()
            if rand > args.rate:
                num_rate_skipped+=1 
                continue
            
            num_mapped += 1

            called_seq, aligned_seq = get_called_and_aligned(read, genome_mask, skip_ins=True)

            # Read is skipped if cigar string contains operations other than CMATCH, CSOFT_CLIP, CINS, CDEL
            if not called_seq:
                log.debug("Cigar string doesn't meet requirements, skipping read")
                continue

            byte_called_seq = byte_encode_seq(called_seq, max_read_length=args.read_length)
            byte_aligned_seq = byte_encode_seq(aligned_seq, max_read_length=args.read_length)
            
            assert byte_called_seq.shape == byte_aligned_seq.shape
            no_call_mask = np.logical_and(byte_called_seq != NUCLEOTIDE2IDX['N'], byte_aligned_seq != NUCLEOTIDE2IDX['N'])
            del_mask = byte_called_seq != NUCLEOTIDE2IDX['D']
            variant_mask = byte_called_seq != NUCLEOTIDE2IDX['V']
            relevant_basecalls_mask = np.logical_and.reduce(np.array([no_call_mask, del_mask, variant_mask]))
            
            mismatch_seq = np.zeros_like(byte_called_seq, dtype=bool)            
            mismatch_seq[relevant_basecalls_mask] = byte_called_seq[relevant_basecalls_mask] != byte_aligned_seq[relevant_basecalls_mask]

            # iterate over all kmers
            n_kmers = len(aligned_seq) - args.kmer_length + 1
            
            for i in range(n_kmers):
                kmer = aligned_seq[i:i + args.kmer_length]
                mask = relevant_basecalls_mask[i:i + args.kmer_length]
                mismatch = mismatch_seq[i:i + args.kmer_length]

                # skip any kmer overlapping a D/N/V masked position
                if np.sum(mask) != args.kmer_length:
                    log.debug(f"skipping kmer {kmer}, {np.sum(mask)} != {args.kmer_length}")
                    continue
                    
                num_kmers_used += 1
                
                kmer_dict[kmer]["all"] = np.add(kmer_dict[kmer]["all"],np.ones(len(kmer),dtype=int))
                kmer_dict[kmer]["mismatch"] = np.add(kmer_dict[kmer]["mismatch"],mismatch)
                if np.sum(mismatch) == 0:
                    kmer_dict[kmer]["correct"] += 1 
                else:
                    kmer_dict[kmer]["incorrect"] += 1

        bam_fh.close()

    log.info(f"num_secondary={num_secondary}")
    log.info(f"num_supplementary={num_supplementary}")
    log.info(f"num_unmapped={num_unmapped}")
    log.info(f"num_mapped={num_mapped}")
    log.info(f"num_below_mapq_thresh={num_below_mapq_thresh}")
    log.info(f"num_rate_skipped={num_rate_skipped}")
    log.info(f"num_mapped={num_mapped}")
    log.info(f"num_kmers={num_kmers_used}")

    log.info("dumping kmer-error")
    interval_error_tsv = os.path.join(args.output_prefix+".kmer-error.tsv")
    with open(interval_error_tsv, "w") as out_tsv:

        for kmer in kmer_vocabulary:
            kmer_all = kmer_dict[kmer]["all"]
            kmer_mismatch = kmer_dict[kmer]["mismatch"] 
            kmer_correct = kmer_dict[kmer]["correct"]
            kmer_incorrect = kmer_dict[kmer]["incorrect"]
            
            num_kmer_observations = np.max(kmer_all)
            
            assert num_kmer_observations == kmer_correct + kmer_incorrect
            kmer_incorrect_rate = None
            if num_kmer_observations > 0:
                kmer_incorrect_rate = kmer_incorrect / num_kmer_observations

            num_kmer_observations = np.max(kmer_all)
            kmer_base_error_rate = (kmer_mismatch / kmer_all)
            kmer_error_rate = None
            kmer_error_rate = float(np.sum(kmer_mismatch) / np.sum(kmer_all))

            print(kmer, args.chr, num_kmer_observations, ','.join([str(ka) for ka in kmer_all]), ','.join([str(km) for km in kmer_mismatch]), ','.join([str(kbep) for kbep in kmer_base_error_rate]), kmer_error_rate, kmer_correct, kmer_incorrect, kmer_incorrect_rate, sep="\t", file=out_tsv)


if __name__ == "__main__":
    main()
