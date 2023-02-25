# AvidityManuscript2023

Bioinformatic tools and scripts to support Arslan et. al., Sequencing by avidity enables high accuracy with low reagent consumption.

Andrew Altomare
Bryan R Lajoie
Edmund Miller
Ryan Kelley
Semyon Kruglyak

## Abstract
We present avidity sequencing - a novel sequencing chemistry that separately optimizes the process of stepping along a DNA template and the process of identifying each nucleotide within the template.  Nucleotide identification uses multivalent nucleotide ligands on dye-labeled cores to form polymerase-polymer nucleotide complexes bound to clonal copies of DNA targets.  These polymer-nucleotide substrates, termed avidites, decrease the required concentration of reporting nucleotides from micromolar to nanomolar, and yield negligible dissociation rates.  We demonstrate the use of avidites as a key component of a sequencing technology that surpasses Q40 accuracy and enables a diversity of applications that include single cell RNA-seq and whole human genome sequencing.  We also show the advantages of this technology in sequencing through long homopolymers.  

biorxiv - https://www.biorxiv.org/content/10.1101/2022.11.03.514117v1
doi: https://doi.org/10.1101/2022.11.03.514117

<tbd publication link>

## Setup/Install
```
mamba env create -f environment.yml
mamba activate AvidityManuscript2023
```

## Whole Genome Sequencing Analysis

Whole genome sequencing analysis
A FASTQ file with the base calls and quality scores was down-sampled to 35X raw coverage (360,320,126 Input reads) and used as an input into Sentieon BWA following by Sentieon DNAscope [44].  Following alignment and variant calling, the variant calls were compared to the NIST genome in a bottle truth set v4.2.1 via the hap.py comparison framework to derive total error counts and F1 scores[45].  The results are computed based on the 3,848,590 SNV and 982,234 indel passing variant calls made by DNAScope.

## Figure-1
not applicable


## Figure-2
not applicable


## Figure-3

Figure 3: Predicted and observed quality scores for a 2x150 bp sequencing run of human genome HG002.  Panel (a) shows read 1 and the panel (b) shows read 2.  Points on the diagonal indicate that predicted scores match observed scores.  The histograms show that the majority of the data points are above Q40, or 1 error in 10,000 bp.

To assess the accuracy of quality scores shown in Fig. 3, the FASTQ files were aligned with BWA to generate BAM files. GATK BaseRecalibrartor was then applied to the BAM, specifying publicly available known sites files to exclude human variant positions. 

The command used is found below:     

```
gatk BaseRecalibrator --preserve-qscores-less-than 0 -R genome.fa -I sample.bam --known-sites HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites dbsnp_144.hg38.vcf.gz -L HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -O sample.table
```

## Figure-4

Figure 4: The mismatch percentage of AVITI, NovaSeq 6000, and NextSeq 2000 reads before and after homopolymers of length 12 or greater.  

A BED file provided by NIST genome-stratifications v3.0, containing 673,650 homopolymers of length greater than 11 was used to define the regions of interest for the homopolymer analysis (GRCh38_SimpleRepeat_homopolymer_gt11_slop5).

s3://giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_SimpleRepeat_homopolymer_gt11_slop5.bed.gz.  

Reads that overlapped these BED intervals (using samtools view -L and adjusting for the slop5) were selected for accuracy analysis.  Reads with any of the following flags set were discarded (secondary, supplementary, unmapped or reads with mapping quality of 0).  Reads were oriented in the 5’ -> 3’ direction, and split into 3 segments, preceding the homopolymer, overlapping the homopolymer, and following the homopolymer.  The mismatch rate for each read-segment was computed, excluding N-calls, softclipped bases and indels.  For example, if a 150 bp read (aligned on the forward strand) contains a homopolymer in positions 100-120, then the first 99 cycles were used to compute the error rate prior to the homopolymer, and the last 30 cycles were used to compute the error rate following the homopolymer.  Reads were discarded if either the sequence preceding or following the homopolymer was less than 5bp in length.   All reads were then stacked into a matrix, according to their positional offset relative to the homopolymer, and error rate per pos-offset was computed.

The average error rate was computed for avidity sequencing runs and for publicly available data from multiple SBS instruments, for comparison.  The differences of mismatch percentages, across all BED intervals, between AVITI™ and NovaSeq were plotted in a histogram and examples showing various percentiles within the distribution were chosen for display via IGV.

Publicly available data sets for NovaSeq were obtained from the Google Brain Public Data repository on Google Cloud [42].  Publicly available NextSeq 2000 data was obtained from Illumina demo data on BaseSpace [43]. 

The *.interval-error.tsv and *.offset-error.tsv files can be found in the following directory:
https://github.com/Elembio/AvidityManuscript2023/tree/main/data/homopolymer-error/GRCh38_SimpleRepeat_homopolymer_gt11_slop5

The command used to generate the homopoylmer-error:
```
THREADS=6 CONCURRENCY=3 RATE=1.0 SLOP=5 MIN_INTERVAL=1 MIN_PRE=10 MIN_POST=10 bash AvidityManuscript2023/bash/run_stack_reads.sh -r <path_to_sample_specific_dir_containing_bam> -b AvidityManuscript2023/data/bed/GRCh38_SimpleRepeat_homopolymer_gt11_slop5/ -g AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -i test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -o test__GRCh38_SimpleRepeat_homopolymer_gt11_slop5 -w <path_to_scratch_dir>
```

The scripts/notebooks supporting the plots and figure generation:
```
jupyter lab --no-browser
<notebooks/compare_read_stack.ipynb>
```

## Figure-5
 
Figure 5: The mismatch rate comparison following homopolymers lengths 4 through 29.  The mismatch percent difference between avidity sequencing and SBS increases with homopolymer length.  The box plot shows median, quartiles, and the whiskers are 1.5*IQR.  

For this analysis, the following GIAB supplied bed files were combined, and duplicates were removed:

s3://giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed.gz
s3://giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_SimpleRepeat_homopolymer_7to11_slop5.bed.gz
s3://giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_SimpleRepeat_homopolymer_gt11_slop5.bed.gz
s3://giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_SimpleRepeat_homopolymer_gt20_slop5.bed.gz

Producing a new bed file representing all homopolymer of size 4 to inf.
GRCh38_SimpleRepeat_homopolymer_4toinf_slop5

Data was analyzed and plotted as described in Figure 4 (https://github.com/Elembio/AvidityManuscript2023/blob/main/README.md#figure-4)


## ExtendedDataFigure 1 
not applicable

## ExtendedDataFigure 2

Extended Data Fig. 2: Percentage of instances that a k-mer contained at least one mismatch compared across 3 instruments.  Panels A, B, and C display 1-mers, 2-mers, and 3-mers, respectively.  The bars are sorted by AVITI contexts from most to least accurate.      

The same run used to generate the recalibrated quality scores was analyzed via custom script for all k-mers of size 1, 2, and 3.  The computation is based on 1% of a 35X genome to ensure adequate sampling of each k-mer.  For example, each 3-mer is sampled at least 850 thousand times with an average of 6.7 million times.  The figure is based on a publicly available run from each platform.  For the instances of each k-mer, the percent mismatching a variant-masked reference was computed.  The same script was applied to a publicly available NovaSeq data set for HG002 and a publicly available NextSeq 2000 data set for HG001 (demo data for HG002 was not available).  We tabulated the number of k-mers in which the percent incorrect was lowest for AVITI among the three platforms compared.   

The command used to generate the kmer-error across the three K sizes [1,2,3]:
```
#kmer-1
THREADS=6 CONCURRENCY=3 RATE=0.1 SAMPLE=HG002 bash AvidityManuscript2023/bash/run_error_by_kmer.sh -r <path_to_sample_specific_dir_containing_bam> -g AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -b AvidityManuscript2023/data/bed/GRCh38_10bp_1000000_random/ -i test__kmer-1 -o test_kmer-1/ -k 1 -w <path_to_scratch_dir>
#kmer-2
THREADS=6 CONCURRENCY=3 RATE=0.1 SAMPLE=HG002 bash AvidityManuscript2023/bash/run_error_by_kmer.sh -r <path_to_sample_specific_dir_containing_bam> -g AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -b AvidityManuscript2023/data/bed/GRCh38_10bp_1000000_random/ -i test__kmer-2 -o test_kmer-2/ -k 2 -w <path_to_scratch_dir>
#kmer-3
THREADS=6 CONCURRENCY=3 RATE=0.1 SAMPLE=HG002 bash AvidityManuscript2023/bash/run_error_by_kmer.sh -r <path_to_sample_specific_dir_containing_bam> -g AvidityManuscript2023/data/genome/Homo_sapiens_assembly38_primaryonly/ -b AvidityManuscript2023/data/bed/GRCh38_10bp_1000000_random/ -i test__kmer-3 -o test_kmer-3/ -k 3 -w <path_to_scratch_dir>
```

The scripts/notebooks supporting the plots and figure generation:
```
jupyter lab --no-browser
<notebooks/compare_kmer_error.ipynb>
```

## ExtendedDataFigure 3

Extended Data Fig. 3: Histogram of pairwise error differences.  Difference was selected as the metric to cancel the effects of human variants from the mismatch percent.

Data was analyzed and plotted as described in Figure 4 (https://github.com/Elembio/AvidityManuscript2023/blob/main/README.md#figure-4)

## ExtendedDataFigure 4

Extended Data Fig. 4: IGV display of homopolymer loci at the 5thth, 50thth, and 95thth percentile of AVITI minus NovaSeq mismatch percent (corresponding to the dashed lines of Extended Data Fig. 3).  The red bar at the top indicates the homopolymer.  Colors within the IGV read stack correspond to mismatches and softclipping.  Only mismatches contribute to the error rate calculation and softclipped bases are ignored.    

Data was analyzed and plotted as described in Figure 4 (https://github.com/Elembio/AvidityManuscript2023/blob/main/README.md#figure-4)

## ExtendedDataFigure 5 

Extended Data Fig. 5: Comparison of read number vs genomic coverage computed via Picard for PCR-free whole genome data.  AVITI most closely matches the 45-degree line due to the low duplicate rate. 

Another common application is human whole genome sequencing.  This application challenges sequencer accuracy to a greater extent than measuring gene expression because the latter requires only accurate alignment while the former depends on nucleotide accuracy to resolve variant calls.  To demonstrate performance for this application, the well characterized human sample HG002 was prepared for sequencing using a Covaris shearing and PCR-free library preparation method and sequenced with 2x150bp reads.  The run generated 1.02 billion passing filter paired-end reads with a duplicate rate of 0.58% (0.11% classified as optical duplicates by Picard[37]).  To underscore the impact of low duplicates, we compared the number of input reads to genomic coverag (Extended Data Fig. 5).

Data was analyzed as described in WGS )https://github.com/Elembio/AvidityManuscript2023#whole-genome-sequencing-analysis)

## ExtendedDataFigure 6

Extended Data Fig. 6: F1 score for SNPs and indels stratified by all GiaB regions with at least 100 variants in the 4.2.1 truth set of sample HG002.

A FASTQ file with the base calls and quality scores was down-sampled to 35X coverage and used as an input into the DNAScope analysis pipeline from Sentieon.  SNP and indel calls achieved F1 scores of 0.995 and 0.996, respectively.  Table 2 shows variant calling performance for SNPs and small indels on the GIAB-HC regions.  Sensitivity, precision, and F1-score are shown.  The performance on SNPs and indels is comparable.  Extended Data Fig. 6 shows the F1 score for SNPs and indels across all GiaB stratifictions with at least 100 variants in the truth set.  

Data was analyzed as described in WGS )https://github.com/Elembio/AvidityManuscript2023#whole-genome-sequencing-analysis)

## Open Jupyter-Lab notebooks
```
jupyter lab --no-browser
```