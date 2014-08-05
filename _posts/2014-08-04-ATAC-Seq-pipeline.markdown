---
layout: post
title:  "ATAC-Seq Analysis Pipeline"
date:  Mon Aug  4 14:59:01 PDT 2014
categories: orolab ngs atac-seq
---

Generally, the pipeline includes the following steps:

* Align the reads to reference genome by [Bowtie][bowtie]
* Make UCSC track
* Call peaks by [MACS][macs] (use macs2 instead of macs14)
* Annotating peaks by [HOMER][homer]

**Note**, You will need two additional scripts. One is called [getFragmentBed.pl][getFragmentBed] which will convert the mapping coordinates of each mate to the biological fragment level; 
the other is called [preShift.pl][preShift] which shifts the reads mapping to plus strand upstream ***n bp*** and reads mapping to the minus strand downstream ***n bp***.

[bowtie]:  http://bowtie-bio.sourceforge.net/index.shtml
[macs]:    http://liulab.dfci.harvard.edu/MACS/
[getFragmentBed]: https://github.com/riverlee/ATAC/blob/master/code/getFragmentBed.pl
[homer]:   http://biowhat.ucsd.edu/homer/ngs/
[preShift]: https://github.com/riverlee/ATAC/blob/master/code/preShift.pl

{% highlight bash %}
#!/bin/bash
###################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
###################################
# ATAC-Seq analysis pipeline
###################################


###################################
# An example of ATAC-Seq analysis
# Suppose your input file is R1.fq and R2.fq
###################################
###################################
# Alignment
###################################
# 1.1 Alignment 
echo 'bowtie alignment' `date`
bowtie --chunkmbs 256 -p 24 -S  -m 1 -X 2000  -t hg19  R1.fq R2.fq my.sam
samtools view -bS  my.sam |samtools sort - my

# 1.2 Sam to bam and only output those mapped
samtools view -b -F my.bam >my_norm_sorted.bam

# 1.3 Delete not sorted bam
echo "delete bam" `date`
rm -rf my.sam

# 1.4 remove duplicates and them remove chrM reads
echo "remove duplicates" `date`
samtools rmdup  my_norm_sorted.bam my_norm_sorted_rmdup.bam
echo "Remove chrM reads" `date`
samtools view -h my_norm_sorted_rmdup.bam |perl -lane 'print $_ if $F[2] ne "chrM"' |samtools view -bS - >my_norm_sorted_rmdup_nochrM.bam
samtools index my_norm_sorted_rmdup_nochrM.bam

# 1.5 Make Fragment Bed 
getFragmentBed.pl my_norm_sorted_rmdup_nochrM.bam my_fragment.bed

###################################
# Make UCSC Tracks
###################################
# 2.1 make bedGraph
echo "make bedgraph" `date`
sort -k1,1 -k2,2n my_fragment.bed |bedItemOverlapCount hg19 -chromSize=hg19.chrom.sizes stdin |sort -k1,1 -k2,2n >my.bedGraph

# 2.2 make bigwig
bedGraphToBigWig my.bedGraph hg19.chrom.sizes my.bw


###################################
# Call peaks
# Before appling macs2, pre-shift the reads
###################################
preShiftBed.pl my_fragment.bed 75 my_fragment_preShift75.bed

test ! -r callpeaks && mkdir callpeaks
cd callpeaks
macs2 callpeak  -t ../my_fragment_preShift75.bed -f BED  -g hs -n my --nomodel --shiftsize 75

###################################
# Peak annotation by Homme
###################################
annotatePeaks.pl my_peaks.narrowPeak hg19 > my_peaks_anno.txt

{% endhighlight %}


