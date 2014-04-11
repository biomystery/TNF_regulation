nascent <- read.csv('../expdata/nascent.csv')
mRNA_predict <- read.csv('./simData/mRNA_all.csv',header=F)
mRNA_expdata <- read.csv('../expdata/mRNA.csv')

############################################################
# plot fig3s1.pdf

pdf(file='fig3s1.pdf', height=3, width=6, onefile=F, family='Helvetica', paper='special', pointsize=10) 

layout(matrix(c(1,2),1,2,byrow=T))
colors <- c('steelblue4','skyblue1','steelblue3')
pchs <- c(18,17,15)
ylims <- c(0,ceiling(max(nascent[,-1]))+5)

# plot 3: lower process in knockouts. 
matplot(mRNA_predict[,1],mRNA_predict[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (AU) ',
        xlim=c(0,120),lwd=2,main='Simulation')

matlines(mRNA_predict[,1],mRNA_predict[,5:6],type='l',col=colors[2:3],
        lty=rep(2,2),xlim=c(0,120),lwd=2)
matlines(mRNA_predict[,1],mRNA_predict[,7:8],type='l',col=colors[2:3],
        lty=rep(3,2),xlim=c(0,120),lwd=2)


# plot 4: Output exp
#install.packages('plotrix')
library(plotrix)
matplot(mRNA_expdata[,1],mRNA_expdata[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (Normalized)',xlim=c(0,120),
        lwd=2, add =F)
title(main='Experiment',col.main = 'red')
legend('bottom',c("wt","mko","tko"),lty=rep(1,3),pch=pchs,col=colors,bty="n",cex=.8)

# ref : https://stat.ethz.ch/pipermail/r-help/2009-April/195287.html
#plotCI(x=rep(mRNA_expdata[,1],3),y=as.vector(mRNA_expdata[,c(2,4,6)]),
#       uiw=as.vector(mRNA_expdata[,c(3,5,7)]),
#       col=colors,add=T)

dev.off()
system('open fig3s1.pdf')

