import argparse
import gzip
import os
import gc
import sys
import logging
import random
import numpy as np
from tqdm import tqdm
import pysam
import subprocess
from itertools import product

def main():
    
    parser = argparse.ArgumentParser(description='combine compute-error-by-kmer tsv split by region/chr')

    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--kmer_length', type=int, default=8, help='length of kmer to compute error for')
    parser.add_argument('-o', '--output-prefix', type=str, required=True, default=None, help='output prefix - workflow name')

    args = parser.parse_args()

    buffer = {}
    buffer["num_kmer_observations"] = 0
    buffer["kmer_all"] = np.zeros(args.kmer_length,dtype=int)
    buffer["kmer_mismatch"] = np.zeros(args.kmer_length,dtype=int)
    buffer["kmer_correct"] = 0
    buffer["kmer_incorrect"] = 0
    
    kmer_error_tsv = os.path.join(args.output_prefix+".kmer-error.tsv")
    with open(kmer_error_tsv, "w") as out_tsv:
        
        last_kmer = None
        for line in sys.stdin:
            line = line.rstrip()
            fields = line.split("\t")
            
            kmer = fields[0]
            chrom = fields[1]
            num_kmer_observations = fields[2]
            kmer_all = fields[3]
            kmer_mismatch = fields[4]
            kmer_correct = fields[7]
            kmer_incorrect = fields[8]

            if last_kmer is None:
                last_kmer = kmer
            
            if kmer != last_kmer:
                # dump buffer / combine
                combined_num_kmer_observations = buffer["num_kmer_observations"]
                combined_kmer_all = buffer["kmer_all"]
                combined_kmer_mismatch = buffer["kmer_mismatch"]
                combined_kmer_correct = buffer["kmer_correct"]
                combined_kmer_incorrect = buffer["kmer_incorrect"]
                
                assert combined_num_kmer_observations == combined_kmer_correct + combined_kmer_incorrect
                combined_kmer_incorrect_rate = None
                if combined_num_kmer_observations > 0:
                    combined_kmer_incorrect_rate = combined_kmer_incorrect / combined_num_kmer_observations

                combined_num_kmer_observations = np.max(combined_kmer_all)
                combined_kmer_base_error_rate = (combined_kmer_mismatch / combined_kmer_all)
                combined_kmer_error_rate = float(np.sum(combined_kmer_mismatch) / np.sum(combined_kmer_all))

                print(last_kmer, "", combined_num_kmer_observations, ','.join([str(ka) for ka in combined_kmer_all]), ','.join([str(km) for km in combined_kmer_mismatch]), ','.join([str(kbep) for kbep in combined_kmer_base_error_rate]), combined_kmer_error_rate, combined_kmer_correct, combined_kmer_incorrect, combined_kmer_incorrect_rate, sep="\t", file=out_tsv)

                buffer = {}
                buffer["num_kmer_observations"] = 0
                buffer["kmer_all"] = np.zeros(args.kmer_length,dtype=int)
                buffer["kmer_mismatch"] = np.zeros(args.kmer_length,dtype=int)
                buffer["kmer_correct"] = 0
                buffer["kmer_incorrect"] = 0

            #print(kmer, chrom, num_kmer_observations, kmer_all, kmer_mismatch, kmer_correct, kmer_incorrect, sep="\t")
            buffer["kmer"] = kmer
            buffer["chrom"] = ""
            buffer["num_kmer_observations"] += int(num_kmer_observations)
            buffer["kmer_all"] = np.add(buffer["kmer_all"],  np.array(kmer_all.split(","), dtype=int))
            buffer["kmer_mismatch"] = np.add(buffer["kmer_mismatch"], np.array(kmer_mismatch.split(","), dtype=int))
            buffer["kmer_correct"] += int(kmer_correct)
            buffer["kmer_incorrect"] += int(kmer_incorrect)

            last_kmer = kmer

        # dump buffer / combine
        combined_num_kmer_observations = buffer["num_kmer_observations"]
        combined_kmer_all = buffer["kmer_all"]
        combined_kmer_mismatch = buffer["kmer_mismatch"]
        combined_kmer_correct = buffer["kmer_correct"]
        combined_kmer_incorrect = buffer["kmer_incorrect"]

        assert combined_num_kmer_observations == combined_kmer_correct + combined_kmer_incorrect
        combined_kmer_incorrect_rate = None
        if combined_num_kmer_observations > 0:
            combined_kmer_incorrect_rate = combined_kmer_incorrect / combined_num_kmer_observations

        combined_num_kmer_observations = np.max(combined_kmer_all)
        combined_kmer_base_error_rate = (combined_kmer_mismatch / combined_kmer_all)
        combined_kmer_error_rate = float(np.sum(combined_kmer_mismatch) / np.sum(combined_kmer_all))

        print(last_kmer, "", combined_num_kmer_observations, ','.join([str(ka) for ka in combined_kmer_all]), ','.join([str(km) for km in combined_kmer_mismatch]), ','.join([str(kbep) for kbep in combined_kmer_base_error_rate]), combined_kmer_error_rate, combined_kmer_correct, combined_kmer_incorrect, combined_kmer_incorrect_rate, sep="\t", file=out_tsv)

if __name__ == "__main__":
    main()