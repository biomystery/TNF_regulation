# read data

proTNF <- read.csv('proTNF.csv')#,header=F)
facTNF <-read.csv('facs.csv')#,header=F)
pro_scale <- max(proTNF[,2:4])#100
fac_scale <- max(facTNF[,2:4])
combinedTNF <- read.csv('combineTNF.csv')
  
  
# original data 
#http://robjhyndman.com/hyndsight/r-graph-with-two-y-axes/  
par(mar=c(5,4,4,5)+0.1)
matplot(proTNF[,1],proTNF[,2:4],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF',
        ylim=c(-10,100))
?matplot
par(new=T)
matplot(facTNF[,1],facTNF[,2:4]/fac_scale*100,col=c('black','purple','cyan3'),
        type='b',pch=rep(2,3),
        lty=rep(2,3),xlim=c(0,120),xlab='',ylab='',
        xaxt='n',yaxt='n',ylim=c(-10,100))
axis(4)
mtext('facTNF(triangle,dashed lines)',side=4,line=3)
par(new=T)
matplot(combinedTNF[,1],combinedTNF[,2:4],col=c('black','purple','cyan3'),
        type='l',lwd=4,
        lty=rep(1,3),xlim=c(0,120),xlab='',ylab='',
        xaxt='n',yaxt='n',ylim=c(-10,100))
mtext('bold lines: average of two dataset'
      ,cex=1,col='grey')
par(mfrow=c(1,1))

dev.copy2pdf(file="proTNFCombined.pdf", width = 8, height = 7)
system('open proTNFCombined.pdf')
