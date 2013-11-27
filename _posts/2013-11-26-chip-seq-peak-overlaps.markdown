---
layout: post
title:  "Chip-Seq Peak overlaps"
date:   Tue Nov 26 16:37:53 PST 2013
categories: orolab ngs chip-seq
---

I used the following two softwares to do the peaks overlapping. Bedtools is easy and fast while ChIPpeakAnno could provide you more detailed information.

* [bedtools][bedtools]: a flexible suite of utilities for comparing genomic features
* [ChIPpeakAnno][ChIPpeakAnno] - An R pacakge.

[bedtools]: https://code.google.com/p/bedtools/ 
[ChIPpeakAnno]: http://www.bioconductor.org/packages/2.13/bioc/html/ChIPpeakAnno.html

<br/>
## bedtools

Suppose we have peaks from P63 (p63.bed) and Foxo3 (fox3.bed), we would like to see which Foxo3 peaks overlap with P63

{% highlight bash%}
# without gap
intersectBed -a fox3.bed -b p63.bed -wa > fox3.ovlp.p63.bed

# with a 2k gap
windowBed -a fox3.bed -b p63.bed -w 2000 -u >fox3.ovlp.2k.p63.bed

{% endhighlight %}

<br/>
## ChIPpeakAnno

Codes work in R environment.

{% highlight r %}
# load library
library("ChIPpeakAnno")

# Read the bed file and store it as GRangedata object
readbed<-function(f){
  peak.data.frame<-read.table(f,header=FALSE,sep="\t",stringsAsFactors=FALSE)[,1:3]
  peak.range<-BED2RangedData(peak.data.frame,header=FALSE)
  peak.range
}

# Read peaks
foxo3<-readbed("fox3.bed")
p63<-readbed("p63.bed")

# Do 2k gap overlaps
o<-findOverlappingPeaks(foxo3,p63,maxgap=2000,select="first",NameOfPeaks1="Foxo3",NameOfPeaks2="P63",annotate=1)

## Write out overlapping peaks
out.dat<-o$OverlappingPeaks[c(2,4,5,7,8,10)]
out.dat$chr<-paste("chr",out.dat$chr,sep="")
write.table(out.dat,file="fox3.ovlp.2k.p63.bed",sep="\t",quote=FALSE,row.names=FALSE)

## Write out the merge peaks.
out.dat<-as.data.frame(o$MergedPeaks)[1:3]
out.dat$space<-paste("chr",out.dat$space,sep="")
write.table(out.dat,file="fox3.p63.2k.merged.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

{% endhighlight %}


