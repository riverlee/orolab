---
layout: post
title:  "Generating matrix for treeview"
date:   Mon Dec  2 11:12:13 PST 2013
categories: orolab ngs chip-seq
---

This script helps you to visualize occupancy around peak centers in a heatmap way. The input are peaks (either peak center file (chr center) or bed format file (chr start end)) and alignment result (either in bam format or bed format). The script is available on my [github][github].

[github]: https://github.com/riverlee/Personal/blob/master/getCenteredSig.pl 

<br/>

## Usage 

{% highlight bash %}
Usage: getCenteredSig.pl -a <alignfile> -c <centerfile> -o <output> -e [extend bp] -w [window size] -n/-non
        -a  alignment file in bed format
        -c  center position file, format is "chr    pos"
        -b  is the center input file in bed format, thus will estimate the center by (start+end)/2, defalut 0, values is 0 or 1
        -o  output file name
        -e  extend length around the center, extend by upstream extend/2 bp and downstream extend/2bp, default 8000bp
        -w  window size to get the signal, default is 10,(1 means get the average coverage for each base around the extended region)
        -n  normalized the signal to RPKM or use the raw count, values is 0 or 1
        -bam whether input alignment file is in bam/sam format, values is 0 or 1,default 0
        -len the average read length, default is 36
        -h  display this help message

{% endhighlight %}
### Example

Suppose you have peaks from P63 named as p63.bed, and you would like to see the JUN occupancy around these p63 peaks in K562 cell line (the alignment of JUN in K562 is named as K562_JUN.bam).


{% highlight bash %}
perl getCenteredSig.pl -c p63.bed -a K562_JUN.bam -b 1 -o JUN_P63_8kb_w20bp.signal.txt -n 1 -bam 1
{% endhighlight %}

## Visualization

There are two ways to visualize the data.

* By R
* Use [Cluster][Cluster] first and then apply [Treeview][Treeview]

[Cluster]: http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm
[Treeview]: http://jtreeview.sourceforge.net/

### By R

{% highlight r %}
# 1) Load data
my.read<-function(f){
  d<-read.table(f,sep="\t",header=FALSE,row.names=1)
  d<-log2(d+1)
  tmp<-unlist(d)
  q99<-quantile(tmp,probs=0.99)
  d[d>=q99]<-q99
  d
}

dat1<-my.read("JUN_P63_8kb_w20bp.signal.txt")
dat1.rowMean<-rowMeans(dat1)
dat1.colMean<-colMeans(dat1)

# 2) Sort data
myorder<-order(dat1.rowMean,decreasing=FALSE)
dat1<-dat1[myorder,]

# 3) Draw heatmap
min.raw<-round(min(dat1),2)
max.raw<-round(max(dat1),2)
cols<-colorRampPalette(c("green", "black", "red"))(100)

z <- seq(min.raw, max.raw, length = length(cols))

png(filename="JUN_P63_8k.png",width=1600,heigh=1200,res=200)

par(mar=c(2,1,1,1))
layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=TRUE),heights=c(1,4,4))

# 3.1 for the color bar
image(z = matrix(z, ncol = 1), col = cols, xaxt = "n", yaxt = "n")
axis(side=1,at=c(0,1),c(min.raw,max.raw),tick=F)

# 3.2 for the heatmap
image(1:ncol(dat1),1:nrow(dat1),t(dat1),xlab="",ylab="",col=cols,axes=FALSE,main="H1esc",cex.main=0.9)
abline(v=200,col='yellow',lwd=1)

# 3.3 for the average signal around center
ymax<-max(dat1.colMean)
plot(dat1.colMean,xlab="",type="l",ylab="Intensity",cex.lab=0.6,cex.axis=0.6,xaxt='n',col='blue',mgp = c(2, 1, 0),lty=1,ylim=c(0,ymax))
axis(side=1,at=c(1,100,200,300,400,500,600,700,800)/2,labels=c("-4kb","-3kb","-2kb","-1kb","0","+1kb","+2kb","+3kb","+4kb"),cex.lab=0.6,cex.axis=0.6,mgp = c(0, 0, 0
abline(v=200,col='grey',lwd=2)
legend('topright',legend=c("H1esc","Helas3","Hepg2","Huvec","K562"),lty=1:5,col=c('blue','red',"purple","brown","orange"),bty='n')
dev.off()

{% endhighlight %}


### Use Cluster and then appy Treeview

Make use when you use Cluster to do clustering, just cluster rows.

