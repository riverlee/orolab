---
layout: post
title:  "Chip-Seq Analysis Pipeline"
date:   Tue Oct 15 09:46:47 PDT 2013
categories: orolab ngs chip-seq
---

Generally, the pipeline includes the following steps:

* Quality control on raw data by [FastQC][fastqc]
* Align the reads to reference genome by [Bowtie][bowtie]
* Make UCSC track
* Call peaks by [MACS][macs]
* Annotating peaks by [HOMER][homer]

[fastqc]:  http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[bowtie]:  http://bowtie-bio.sourceforge.net/index.shtml
[macs]:    http://liulab.dfci.harvard.edu/MACS/
[homer]:   http://biowhat.ucsd.edu/homer/ngs/

{% highlight bash %}
#!/bin/bash
###################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
###################################
# Chip-Seq analysis pipeline
###################################

###################################
# An example of Chip-Seq analysis
# Suppose your input file is my.fastq
###################################

###################################
# Raw data quality control
###################################
fastqc -o . --nogroup my.fastq
# In this case, we suppose there is no
# further steps to deals with raw data.
# We go ahead alignment with the raw fastq file.

###################################
# Alignment
###################################
# 1.1 Alignment 
 echo 'bowtie alignment' `date`
 bowtie --chunkmbs 256 -p 24 -S -a -m 1 --best --strata -t mm9 my.fastq my.sam 

# 1.2 Sam to bam and only output those mapped
echo "sam to bam (only output mapped)" `date`
samtools view -bS  -F 4 my.sam >my_norm.bam

# 1.3 sort bam
echo "sort bam" `date`
samtools sort my_norm.bam my_norm_sorted

# 1.4 Delete not sorted bam
echo "delete bam" `date`
rm -rf my_norm.bam

# 1.5 remove duplicates
echo "remove duplicates" `date`
samtools rmdup -s my_norm_sorted.bam my_norm_sorted_rmdup.bam

###################################
# Make UCSC Tracks
###################################
# 2.1 make bedGraph
echo "make bedgraph" `date`
genomeCoverageBed -bg -ibam my_norm_sorted_rmdup.bam -split -g mm9.chrom.sizes >my_norm_sorted_rmdup.bedGraph

# 2.2 make bigwig
bedGraphToBigWig my_norm_sorted_rmdup.bedGraph mm9.chrom.sizes my_norm_sorted_rmdup.bw

###################################
# Call peaks
###################################
test ! -r callpeaks && mkdir callpeaks
cd callpeaks
macs14 -t ../my_norm_sorted_rmdup.bam -f BAM -g mm -n my

###################################
# Peak annotation by Homme
###################################
annotatePeaks.pl my_peaks.bed mm9 > my_peaks_anno.txt

{% endhighlight %}

