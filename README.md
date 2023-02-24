# AvidityManuscript2023

Bioinformatic tools and scripts to support Arslan et. al., Sequencing by avidity enables high accuracy with low reagent consumption.
<tbd paper link>

## Setup/Install
```
<download and install conda>
<install samtools, install parallel>
..
conda create -n AvidityManuscript2023 python=3.9
conda activate AvidityManuscript2023
pip3 install -r requirements.txt
```

## Figure-1
```
not applicable
```

## Figure-2
```
not applicable
```

## Figure-3
```
Figure 3: Predicted and observed quality scores for a 2x150 bp sequencing run of human genome HG002.  The left panel shows read 1 and the right panel shows read 2.  Points on the diagonal indicate that predicted scores match observed scores.  The histograms show that the majority of the data points are above Q40, or 1 error in 10,000 bp.

BBS-0174, processed through bwa using `Homo_sapiens_assembly38`.  GATK BaseRecalibrator

`gatk BaseRecalibrator --preserve-qscores-less-than 0 -R genome.fa -I sample.bam --known-sites HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites dbsnp_144.hg38.vcf.gz         -L HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -O sample.table`

## Figure-4
## Figure-5
## Figure-6

## Open Jupyter-Lab notebooks
```
jupyter lab --no-browser
```

## Compute error by K-mer

# small bam
THREADS=10 CONCURRENCY=10 RATE=0.1 SAMPLE=HG002 bash ~/git/AvidityManuscript2023/bash/run_error_by_kmer.sh -r ~/git/AvidityManuscript2023/data/bam/ac9abe80f5044b52bf42c2a40cbd9b6a/ -g ~/git/AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test -o test_output/ -k 3 -b ~/git/AvidityManuscript2023/data/bed/GRCh38_10bp_1000000_random/ -w ~/scratch/kmer/

THREADS=12 CONCURRENCY=12 RATE=1.0 SLOP=0 MIN_INTERVAL=1 MIN_PRE=10 MIN_POST=10 bash ~/git/AvidityManuscript2023/bash/run_stack_reads.sh -r ~/git/AvidityManuscript2023/data/bam/ac9abe80f5044b52bf42c2a40cbd9b6a/ -b ~/git/AvidityManuscript2023/data/bed/GRCh38_10bp_1000000_random/ -g ~/git/AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test__GRCh38_10bp_1000000_random -o test__GRCh38_10bp_1000000_random -w ~/scratch/read_stack__random/

THREADS=12 CONCURRENCY=12 RATE=1.0 SLOP=0 MIN_INTERVAL=1 MIN_PRE=10 MIN_POST=10 bash ~/git/AvidityManuscript2023/bash/run_stack_reads.sh -r ~/git/AvidityManuscript2023/data/bam/ac9abe80f5044b52bf42c2a40cbd9b6a/ -b ~/git/AvidityManuscript2023/data/bed/GRCh38_SimpleRepeat_homopolymer_gt11_slop5/ -g ~/git/AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -o test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -w ~/scratch/read_stack__homopolymer/

# BBS-0174 bam
THREADS=12 CONCURRENCY=12 RATE=0.1 SAMPLE=HG002 bash ~/git/AvidityManuscript2023/bash/run_error_by_kmer.sh -r ~/git/AvidityManuscript2023/data/bam/20220601_PLT-03_BBS-0174-OBPA__FQD-2x150-35x/ -g ~/git/AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test -o test_output/ -k 3 -b ~/git/AvidityManuscript2023/data/bed/GRCh38_10bp_1000000_random/ -w ~/scratch/kmer/

THREADS=12 CONCURRENCY=12 RATE=1.0 SLOP=5 MIN_INTERVAL=1 MIN_PRE=10 MIN_POST=10 bash ~/git/AvidityManuscript2023/bash/run_stack_reads.sh -r ~/git/AvidityManuscript2023/data/bam/20220601_PLT-03_BBS-0174-OBPA__FQD-2x150-35x/ -b ~/git/AvidityManuscript2023/data/bed/GRCh38_SimpleRepeat_homopolymer_gt11_slop5/ -g ~/git/AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -o test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -w ~/scratch/read_stack__homopolymer/