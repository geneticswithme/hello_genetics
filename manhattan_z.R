#主要是gettinggeneticsdone(blogspot)的manhattanplot的Rscript，zyy对其稍加修改，可联系yyzhang@cau.edu.cn
#改动后，CHR可以为连续或不连续的数值或字符。
#使用方法
#对plink的关联分析结果做manhattan图
#plink --file input --assoc 生成的plink.assoc
#打开R，
#df<-read.table("plink.assoc")
#source("manhattan_z.R")
#manhattan(df,col=c(1,2),limitchromosomes=1:23,genomewideline=8,suggestiveline=9)



manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    
    d=d[d$CHR!=0, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    chrnames=unique(d$CHR)
    numchroms=length(unique(d$CHR))
	colors <- rep(colors,numchroms)[1:numchroms]
    if (ymax=="max") ymax<-ceiling(max(d$logp))
    if (ymax<8) ymax<-8
    
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        for (i in 1:numchroms) {
          if (i==1) {
    			d[d$CHR==chrnames[i], ]$pos=d[d$CHR==chrnames[i], ]$BP
    		} else {
    			lastbase=lastbase+tail(subset(d,CHR==chrnames[i-1])$BP, 1)
    			d[d$CHR==chrnames[i], ]$pos=d[d$CHR==chrnames[i], ]$BP+lastbase
    		}
    		ticks=c(ticks, d[d$CHR==chrnames[i], ]$pos[floor(length(d[d$CHR==chrnames[i], ]$pos)/2)+1])
    	}
    }
    
    if (numchroms==1) {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }	else {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=chrnames, ...)
        icol=1
        for (i in 1:numchroms) {
            with(d[d$CHR==chrnames[i], ],points(pos, logp, col=colors[icol], ...))
            icol=icol+1
    	}
    }
    
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="green3", ...)) 
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
}
