# read data

proTNF_exp <- read.csv('combineTNF.csv')
proTNF_sim <- read.csv('proTNF_all.csv',header=F)
proTNF_sim_tl <- read.csv('proTNF_all_tl.csv',header=F)
proTNF_sim_sec <- read.csv('proTNF_all_sec.csv',header=F)
proTNF_sim_both <- read.csv('proTNF_all_both.csv',header=F)
colnames(proTNF_sim)<-colnames(proTNF_exp)
colnames(proTNF_sim_tl)<-colnames(proTNF_exp)
colnames(proTNF_sim_sec)<-colnames(proTNF_exp)
colnames(proTNF_sim_both)<-colnames(proTNF_exp)

secTNF_exp <- read.csv('TNF_secrection.csv')
secTNF_sim <- read.csv('sec_all.csv')
secTNF_sim_tl <- read.csv('sec_all_tl.csv')
secTNF_sim_sec <- read.csv('sec_all_sec.csv')
secTNF_sim_both <- read.csv('sec_all_both.csv')
colnames(secTNF_sim)<-colnames(proTNF_exp)
colnames(secTNF_sim_tl)<-colnames(proTNF_exp)
colnames(secTNF_sim_sec)<-colnames(proTNF_exp)
colnames(secTNF_sim_both)<-colnames(proTNF_exp)

mRNA_exp <- read.csv('mRNA.csv')

########
# plot
########
pdf(file='fig4.pdf', height=11, width=8.5, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12) 

layout(matrix(1:10,5,2,byrow=T))

# no regulation
matplot(proTNF_sim[,1],proTNF_sim[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (simulation)',xlim=c(0,120),ylim=c(0,100),
        main = 'No tl and sec regulation')

legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=c('black','purple','cyan3'),bty="n")
matplot(secTNF_sim[,1],secTNF_sim[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (simulation)',xlim=c(0,120),ylim=c(0,3.5),
        main = 'No tl and sec regulation')

# tl only

matplot(proTNF_sim_tl[,1],proTNF_sim_tl[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (simulation)',xlim=c(0,120),ylim=c(0,100),
        main = 'Tl regulation only')

matplot(secTNF_sim_tl[,1],secTNF_sim_tl[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (simulation)',xlim=c(0,120),ylim=c(0,3.5),
        main = 'Tl regulation only')

# sec only
matplot(proTNF_sim_sec[,1],proTNF_sim_sec[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (simulation)',xlim=c(0,120),ylim=c(0,100),
        main = 'Sec regulation only')

matplot(secTNF_sim_sec[,1],secTNF_sim_sec[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (simulation)',xlim=c(0,120),ylim=c(0,3.5),
        main = 'Sec regulation only')
# both

matplot(proTNF_sim_both[,1],proTNF_sim_both[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (simulation)',xlim=c(0,120),ylim=c(0,100),
        main = 'Tl and sec regulation')

matplot(secTNF_sim_both[,1],secTNF_sim_both[,2:4],type='l',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (simulation)',xlim=c(0,120),ylim=c(0,3.5),
        main = 'Tl and sec regulation')

# exp
matplot(proTNF_exp[,1],proTNF_exp[,2:4],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (experiment)')
title(main = "Experimental data",col.main='red')
matplot(secTNF_exp[,1],secTNF_exp[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (experiment)',xlim=c(0,120))
title(main = "Experimental data",col.main='red')
# ending 


dev.off()
system('open fig4.pdf')

