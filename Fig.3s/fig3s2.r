
score <- read.csv('./simData/fig3s2.csv',header=F)
bestfit <- read.csv('./simData/best_fit.csv',header=F)
expdata <- read.csv('./simData/exp_fit.csv',header=F)
scoreMat <- do.call(rbind,score)
rownames(scoreMat) <- seq(0.1,3,by=0.1)
colnames(scoreMat) <- seq(1,11,by=0.5)
#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }


############################################################
# plot fig3s1.pdf

pdf(file='fig3s2_1.pdf', height=3, width=3, onefile=F, family='Helvetica', paper='special', pointsize=10)

# best fit
colors <- c('steelblue4','skyblue1','steelblue3')
matplot(bestfit[,1],bestfit[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (AU) ',
        xlim=c(0,120),lwd=2,main='Best fit')
pchs <- c(18,17,15)
matpoints(expdata[,1],expdata[,c(2,4,6)],type = "p",col=colors,pch=pchs,
        xlim=c(0,120),cex=1.5)
legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=colors,bty="n")
#plotCI(x=rep(expdata[,1],3),y=as.vector(expdata[,c(2,4,6)]),
#       uiw=as.vector(expdata[,c(3,5,7)]),
#       col=colors,add=T)

dev.off()
system('open fig3s2_1.pdf')

pdf(file='fig3s2_2.pdf', height=4, width=4, onefile=F, family='Helvetica', paper='special', pointsize=10)

# heatmap
breaks=seq(min(score), min(score)+0.5, length.out=21) #41 values

mycol <- colorRampPalette(brewer.pal(11, "RdBu"))(length(breaks)-1)

heatmap.2(scoreMat,
          Rowv=NA, 
          Colv=NA, 
          col=mycol,
          scale='none',
          dendrogram = 'none',
          trace= 'none',
          breaks=breaks,
          xlab='mko_fold',
          ylab='tko_fold',
          density.info="none")
#          lmat=rbind(c(0, 3,4), c(2,1,0)),
#          lwid=c(1.5,4,2))



dev.off()
system('open fig3s2_2.pdf')

