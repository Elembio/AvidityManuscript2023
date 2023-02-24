# AvidityManuscript2023

Bioinformatic tools and scripts to support Arslan et. al., Sequencing by avidity enables high accuracy with low reagent consumption.

## Abstract
We present avidity sequencing - a novel sequencing chemistry that separately optimizes the process of stepping along a DNA template and the process of identifying each nucleotide within the template.  Nucleotide identification uses multivalent nucleotide ligands on dye-labeled cores to form polymerase-polymer nucleotide complexes bound to clonal copies of DNA targets.  These polymer-nucleotide substrates, termed avidites, decrease the required concentration of reporting nucleotides from micromolar to nanomolar, and yield negligible dissociation rates.  We demonstrate the use of avidites as a key component of a sequencing technology that surpasses Q40 accuracy and enables a diversity of applications that include single cell RNA-seq and whole human genome sequencing.  We also show the advantages of this technology in sequencing through long homopolymers.  

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
not applicable


## Figure-2
not applicable


## Figure-3

Figure 3: Predicted and observed quality scores for a 2x150 bp sequencing run of human genome HG002.  The left panel shows read 1 and the right panel shows read 2.  Points on the diagonal indicate that predicted scores match observed scores.  The histograms show that the majority of the data points are above Q40, or 1 error in 10,000 bp.

To assess the accuracy of quality scores shown in Figure 3, the FASTQ files were aligned with BWA to generate BAM files. GATK BaseRecalibrartor was then applied to the BAM, specifying publicly available known sites files to exclude human variant positions. 
The command used is found below:     

```
gatk BaseRecalibrator --preserve-qscores-less-than 0 -R genome.fa -I sample.bam --known-sites HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites dbsnp_144.hg38.vcf.gz -L HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -O sample.table
```

## Figure-4

Figure 4: The mismatch percentage of AVITI and NovaSeq reads before and after homopolymers of length 12 or greater.

A BED file provided by NIST genome-stratifications v3.0, containing 673,650 homopolymers of length greater than 11 was used to define the regions of interest for the homopolymer analysis (s3://giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_SimpleRepeat_homopolymer_gt11_slop5.bed.gz).  Reads that overlapped these BED intervals (using samtools view -L and adjusting for the slop5) were selected for accuracy analysis.  Reads with any of the following flags set were discarded (secondary, supplementary, unmapped or reads with mapping quality of 0).  Reads were oriented in the 5’ -> 3’ direction, and split into 3 segments, preceding the homopolymer, overlapping the homopolymer, and following the homopolymer.  The mismatch rate for each read-segment was computed, excluding N-calls, softclipped bases and indels.  For example, if a 150 bp read (aligned on the forward strand) contains a homopolymer in positions 100-120, then the first 99 cycles were used to compute the error rate prior to the homopolymer, and the last 30 cycles were used to compute the error rate following the homopolymer.  Reads were discarded if either the sequence preceding or following the homopolymer was less than 5bp in length.   All reads were then stacked into a matrix, according to their positional offset relative to the homopolymer, and error rate per pos-offset was computed.

The average error rate was computed for avidity sequencing runs and for publicly available data from multiple SBS instruments, for comparison.  The differences of mismatch percentages, across all BED intervals, between AVITI™ and NovaSeq were plotted in a histogram and examples showing various percentiles within the distribution were chosen for display via IGV.

Publicly available data sets for HiSeq and NovaSeq were obtained from the Google Brain Public Data repository on Google Cloud and from the PrecisionFDA Truth Challenge 2 data repository[41, 42]. 

The *.interval-error.tsv and *.offset-error.tsv files can be found in the below directory:
https://github.com/Elembio/AvidityManuscript2023/tree/main/data/homopolymer-error/GRCh38_SimpleRepeat_homopolymer_gt11_slop5

The command used to generate the homopoylmer-error can be found below.

```
THREADS=12 CONCURRENCY=12 RATE=1.0 SLOP=5 MIN_INTERVAL=1 MIN_PRE=10 MIN_POST=10 bash AvidityManuscript2023/bash/run_stack_reads.sh -r <path_to_sample_specific_dir_containing_bam> -b AvidityManuscript2023/data/bed/GRCh38_SimpleRepeat_homopolymer_gt11_slop5/ -g AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -o test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -w <path_to_scratch_dir>
```

The scripts/notebooks supporting the plots and figure generation can be found below:
```
jupyter lab --no-browser
<notebooks/compare_read_stack.ipynb>
```

## Figure-5
 
Figure 5: 


## ExtendedDataFigure 1 
not applicable

## ExtendedDataFigure 2
not applicable

## ExtendedDataFigure 3
not applicable

## ExtendedDataFigure 4
not applicable

## ExtendedDataFigure 5 
not applicable

## ExtendedDataFigure 6
not applicable

## Open Jupyter-Lab notebooks
```
jupyter lab --no-browser
```