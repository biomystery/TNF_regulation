# read data

proTNF_exp <- read.csv('combineTNF.csv')
proTNF_sim <- read.csv('proTNF_all.csv',header=F)
colnames(proTNF_sim)<-colnames(proTNF_exp)
head(proTNF_sim)
secTNF_exp <- read.csv('TNF_secrection.csv')
secTNF_exp
secTNF_sim <- read.csv('sec_all.csv')
colnames(secTNF_sim)<-colnames(proTNF_exp)
head(secTNF_sim)  

mRNA_exp <- read.csv('mRNA.csv')

# plot

pdf(file='fig4.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12) 

layout(matrix(c(1,0,2,3,4,5),3,2,byrow=T))

matplot(mRNA_exp[,1],mRNA_exp[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (experiment)',xlim=c(0,120))


matplot(proTNF_sim[,1],proTNF_sim[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (simulation)',xlim=c(0,120))

legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=c('black','purple','cyan3'),bty="n")
matplot(secTNF_sim[,1],secTNF_sim[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (simulation)',xlim=c(0,120))

matplot(proTNF_exp[,1],proTNF_exp[,2:4],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (experiment)')

matplot(secTNF_exp[,1],secTNF_exp[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (experiment)',xlim=c(0,120))
?legend
# ending 


dev.off()
system('open fig4.pdf')

