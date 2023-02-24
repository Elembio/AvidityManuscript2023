# PrimaryAnalysisBamStats

## Setup

conda create -n AvidityManuscript2023 python=3.9
conda activate AvidityManuscript2023
pip3 install -r requirements.txt


## Compute error by K-mer

THREADS=10 CONCURRENCY=10 RATE=0.1 SAMPLE=HG002 bash ~/git/AvidityManuscript2023/bash/run_error_by_kmer.sh -r ~/git/AvidityManuscript2023/data/bam/ac9abe80f5044b52bf42c2a40cbd9b6a/ -g ~/git/AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test -o test_output/ -k 3 -b ~/git/AvidityManuscript2023/data/bed/GRCh38_10bp_1000000_random/ -w ~/scratch/kmer/

THREADS=12 CONCURRENCY=12 RATE=1.0 SLOP=0 MIN_INTERVAL=1 MIN_PRE=10 MIN_POST=10 bash ~/git/AvidityManuscript2023/bash/run_stack_reads.sh -r ~/git/AvidityManuscript2023/data/bam/ac9abe80f5044b52bf42c2a40cbd9b6a/ -b ~/git/AvidityManuscript2023/data/bed/GRCh38_10bp_1000000_random/ -g ~/git/AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test__GRCh38_10bp_1000000_random -o test__GRCh38_10bp_1000000_random -w ~/scratch/read_stack__random/

THREADS=12 CONCURRENCY=12 RATE=1.0 SLOP=0 MIN_INTERVAL=1 MIN_PRE=10 MIN_POST=10 bash ~/git/AvidityManuscript2023/bash/run_stack_reads.sh -r ~/git/AvidityManuscript2023/data/bam/ac9abe80f5044b52bf42c2a40cbd9b6a/ -b ~/git/AvidityManuscript2023/data/bed/GRCh38_SimpleRepeat_homopolymer_gt11_slop5/ -g ~/git/AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -o test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -w ~/scratch/read_stack__homopolymer/

