import sys
import gzip
import argparse
import numpy as np
from s3path import S3Path
from pathlib import Path

if sys.platform != 'win32':
    from pybedtools import BedTool

from .clocal_path import ClocalPath


SAMPLEDATA = {
    'HG001': {
        'GENOME': 'hg38',
        'BEDFILES': [
            's3://giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed'
        ],
        'VCFFILES': [
            's3://giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
        ]
    },
    'HG002': {
        'GENOME': 'hg38',
        'BEDFILES': [
            's3://giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed'
        ],
        'VCFFILES': [
            's3://giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
        ]
    },
}

GENOME2DATA = {
    'hg38': {
        'LENGTH': 2875001522,
        'CONTIG2POSITION': {
            'chr1': 0, 'chr2': 248956422, 'chr3': 491149951, 'chr4': 689445510, 'chr5': 879660065,
            'chr6': 1061198324, 'chr7': 1232004303, 'chr8': 1391350276, 'chr9': 1536488912,
            'chr10': 1674883629, 'chr11': 1808681051, 'chr12': 1943767673, 'chr13': 2077042982,
            'chr14': 2191407310, 'chr15': 2298451028, 'chr16': 2400442217, 'chr17': 2490780562,
            'chr18': 2574038003, 'chr19': 2654411288, 'chr20': 2713028904, 'chr21': 2777473071,
            'chr22': 2824183054
        }
    } 
}

class HumanExclusionIntervals:

    def __init__(self, contig2position, exclusion_intervals):
        self.contig2position = contig2position
        self.exclusion_intervals = exclusion_intervals

    def __iter__(self):
        return iter(self.exclusion_intervals)

    def __next__(self):
        return next(self.exclusion_intervals)

    @classmethod
    def from_vcf_gz(cls, filename, contig2position):
        exclusion_intervals = []
        filepath = ClocalPath.construct_path(filename)
        vcf_chr_lens = np.zeros(22)

        with filepath.open(mode='rb') as vcf_gz_handle:
            with gzip.open(vcf_gz_handle, 'rt') as vcf_handle:
                for line in vcf_handle:
                    if line.startswith('#'):
                        # get contig lengths for chr1-22
                        if line.startswith('##contig=<'):
                            split_line = line[13:].split(',')
                            contig_name = split_line[0]
                            if contig_name not in contig2position:
                                continue
                            chr_idx = int(contig_name[3:]) - 1
                            assert split_line[1].startswith('length='), 'cant find chromosome length in vcf file'
                            vcf_chr_lens[chr_idx] = int(split_line[1][7:])
                        else:
                            continue
                    else:
                        # validate chr lengths
                        cum_len = 0
                        for chr_num in range(1, 23):
                            chr_name = f'chr{chr_num}'
                            assert contig2position[chr_name] == cum_len, 'mismatch between expected and vcf contig lengths'
                            cum_len += vcf_chr_lens[chr_num-1]

                        row = line.split('\t')
                        if len(row) != 10:
                            raise RuntimeError('unable to split vcf row')
                        if row[0] not in contig2position:
                            continue
                        # skip any variant that includes indels
                        alt_lens = [len(alt) for alt in row[4].split(',')]
                        if len(row[3]) not in alt_lens:
                            continue

                        chromosome_offset = contig2position[row[0]]
                        ref_location = int(row[1])
                        variant_length = len(row[3])
                        start = chromosome_offset + ref_location - 1 # -1 for 0-based position
                        exclusion_intervals.append( (start, start + variant_length) )

        return cls(contig2position, exclusion_intervals)

    @classmethod
    def from_bed(cls, filename, contig2position):
        filepath = ClocalPath.construct_path(filename)
        exclusion_intervals = []
        interval_start = 0
        prev_contig = None
        with filepath.open(mode='r') as bed_handle:
            bed_file = BedTool(bed_handle)
            for feature in bed_file:
                contig_name = feature.chrom

                # ensure bed file contig ordering matches contig2position
                if prev_contig is None:
                    prev_contig = contig_name
                if contig_name != prev_contig:
                    prev_chr_num = int(prev_contig[3:])
                    chr_num = int(contig_name[3:])
                    assert chr_num == prev_chr_num + 1
                    prev_contig = contig_name

                if contig_name not in contig2position:
                    continue

                # get exclusion interval
                contig_offset = contig2position[contig_name]
                interval_end = contig_offset + feature.start
                assert interval_end > interval_start, 'bed file intervals not sorted'
                exclude_interval = (interval_start, interval_end)
                exclusion_intervals.append(exclude_interval)
                interval_start = contig_offset + feature.stop
                assert interval_start >= interval_end, 'bed file intervals not sorted'
        
        assert len(exclusion_intervals) > 0, 'unable to parse exclusion intervals from bed file'
        
        return cls(contig2position, exclusion_intervals)


class GenomeMask:
    """Bitmask over genome indicating which regions to omit
    """
    def __init__(self, mask, genome):
        self.genome = genome
        self._genome_data = GENOME2DATA[genome]
        self._mask = mask
        self.genome_len = self._genome_data['LENGTH']
        self.contig2position = self._genome_data['CONTIG2POSITION']


    @classmethod
    def from_sample_name(cls, sample_name):
        if sample_name not in SAMPLEDATA or SAMPLEDATA[sample_name]['GENOME'] not in GENOME2DATA:
            raise RuntimeError(f'Cannot generate genome mask for {sample_name}')

        sample_data = SAMPLEDATA[sample_name]
        genome = sample_data['GENOME']
        genome_data = GENOME2DATA[genome]
        genome_length = genome_data['LENGTH']
        contig2position = genome_data['CONTIG2POSITION']
        mask_length = (genome_length + 63) >> 6
        mask = np.zeros(mask_length, dtype=np.uint64)
        obj = cls(mask, genome)

        for filename in sample_data['VCFFILES']:
            print(f'Reading VCF {ClocalPath.construct_path(filename).name}...', end='', flush=True)
            vcf = HumanExclusionIntervals.from_vcf_gz(filename, contig2position)
            print('Done.')
            print('Updating genome mask...', end='', flush=True)
            obj._update_mask(vcf)
            print('Done.')
        for filename in sample_data['BEDFILES']:
            print(f'Reading BED {ClocalPath.construct_path(filename).name}...', end='', flush=True)
            bed = HumanExclusionIntervals.from_bed(filename, contig2position)
            print('Done.')
            print(f'Updating genome mask...', end='', flush=True)
            obj._update_mask(bed)
            print('Done.')
        return obj

    def _update_mask(self, intervals):
        for interval in intervals:
            assert interval[1] < self.genome_len and interval[0] >= 0 and interval[1] > interval[0], f'invalid interval: {interval}'
            first_chunk_idx = interval[0] >> 6
            first_chunk_bit = interval[0] - (first_chunk_idx << 6)
            last_chunk_idx = interval[1] >> 6
            last_chunk_bit = interval[1] - (last_chunk_idx << 6)
            curr_chunk_idx = first_chunk_idx
            while curr_chunk_idx <= last_chunk_idx:
                curr_mask_chunk = self._mask[curr_chunk_idx]
                if curr_chunk_idx == first_chunk_idx:
                    start_bit = first_chunk_bit
                else:
                    start_bit = 0
                if curr_chunk_idx == last_chunk_idx:
                    stop_bit = last_chunk_bit
                else:
                    stop_bit = 64
                set_bits = (1 << (stop_bit - start_bit)) - 1
                disj_operand = np.uint64(set_bits << start_bit) # little endian ordering
                curr_mask_chunk = curr_mask_chunk | disj_operand
                self._mask[curr_chunk_idx] = curr_mask_chunk
                curr_chunk_idx += 1


    def slice(self, contig, start, length):
        linear_position = self.contig2position[contig] + start
        curr_idx = linear_position >> 6
        curr_bit = linear_position - (curr_idx << 6)
        mask_slice = []
        mask_chunk = self._mask[curr_idx]
        while len(mask_slice) < length:
            mask_slice.append((np.int64(mask_chunk) >> curr_bit) & 1)
            curr_bit += 1
            if curr_bit == 64:
                curr_bit = 0
                curr_idx += 1
                mask_chunk = self._mask[curr_idx]
        return mask_slice

    def at(self, idx):
        curr_idx = idx >> 6
        curr_bit = idx - (curr_idx << 6)
        mask_chunk = self._mask[curr_idx]
        return (int(mask_chunk) >> int(curr_bit)) & 1
