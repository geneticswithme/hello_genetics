# GWAS Q-Q Plot
#exdat为一系列SNP的P值文件。
#第一列为SNP名称，第二列为染色体号，第三列是物理位置，第四列为Pvalue。文件没有header。
pvals <- read.table("exdat", header = F)
observed_p <- sort(pvals[,4]) #抽取计算的SNP位点的Pvalue
log_observed_p <- -(log10(observed_p)) #计算得到的P值取负log10对数。

# 构建期望P值的-log10对数。log10括号里分母+1，使得结果不出现0值。
log_expected_p <- -(log10(1:length(observed_p)/(length(observed_p)+1)))

#作图，保存为pdf格式。更多设置参见?pdf。
pdf("qq_plot.pdf")
plot ( c(0,7), c(0,7),col = 'red', lwd =3, type = 'l', xlab = "expected (-logP)",
	ylab = "Observed (-logP)", xlim = c(0,7), ylim = c(0,7), las = 1, xaxs = "i", bty = "l")
points( log_expected_p, log_observed_p, pch =20, cex = 0.5)
dev.off()


#计算genomic control的Inflation λ。
res_p <- observed_p[which(!is.na(observed_p))]
library (snpStats)
pdf (qq_plot2.pdf)
qq.chisq (-2 * log (res_p), df=2, pvals = TRUE, overdisp = FALSE, thin = c(0.8, 1000))# full model test，df=2
#屏幕打印出lambda的结果
dev.off()
