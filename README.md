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
not applicable


## Figure-2
not applicable


## Figure-3

Figure 3: Predicted and observed quality scores for a 2x150 bp sequencing run of human genome HG002.  The left panel shows read 1 and the right panel shows read 2.  Points on the diagonal indicate that predicted scores match observed scores.  The histograms show that the majority of the data points are above Q40, or 1 error in 10,000 bp.

BBS-0174, processed through bwa using `Homo_sapiens_assembly38`.  GATK BaseRecalibrator

```
gatk BaseRecalibrator --preserve-qscores-less-than 0 -R genome.fa -I sample.bam --known-sites HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites dbsnp_144.hg38.vcf.gz         -L HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -O sample.table
```

## Figure-4
not applicable

## Figure-5
not applicable

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