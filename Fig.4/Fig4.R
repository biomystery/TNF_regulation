# read data

proTNF_exp <- read.csv('../expdata/combineTNF.csv')
proTNF_sim <- read.csv('./simData/proTNF_all.csv',header=F)
proTNF_sim_tl <- read.csv('./simData/proTNF_all_tl.csv',header=F)
proTNF_sim_sec <- read.csv('./simData/proTNF_all_sec.csv',header=F)
proTNF_sim_both <- read.csv('./simData/proTNF_all_both.csv',header=F)
colnames(proTNF_sim)<-colnames(proTNF_exp)
colnames(proTNF_sim_tl)<-colnames(proTNF_exp)
colnames(proTNF_sim_sec)<-colnames(proTNF_exp)
colnames(proTNF_sim_both)<-colnames(proTNF_exp)

secTNF_exp <- read.csv('../expdata/TNF_secrection.csv')
secTNF_sim <- read.csv('./simData/sec_all.csv')
secTNF_sim_tl <- read.csv('./simData/sec_all_tl.csv')
secTNF_sim_sec <- read.csv('./simData/sec_all_sec.csv')
secTNF_sim_both <- read.csv('./simData/sec_all_both.csv')
colnames(secTNF_sim)<-colnames(proTNF_exp)
colnames(secTNF_sim_tl)<-colnames(proTNF_exp)
colnames(secTNF_sim_sec)<-colnames(proTNF_exp)
colnames(secTNF_sim_both)<-colnames(proTNF_exp)

mRNA_exp <- read.csv('../expdata/mRNA.csv')

########
# plot
########
pdf(file='fig4.pdf', height=4, width=8, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12) 
colors <- c('steelblue4','skyblue1','steelblue3')
pchs <- c(18,17,15)

layout(matrix(1:6,2,3,byrow=T))

# no regulation
matplot(proTNF_sim[,1],proTNF_sim[,2:4],type='l',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (sim)',xlim=c(0,120),ylim=c(0,100),
        main = 'No tl regulation',lwd=2)

matplot(proTNF_sim_tl[,1],proTNF_sim_tl[,2:4],type='l',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (sim)',xlim=c(0,120),ylim=c(0,100),
        main = 'Tl regulation',lwd=2)



matplot(proTNF_exp[,1],proTNF_exp[,2:4],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='proTNF (exp)',lwd=2)
title(main = "Experimental data",col.main='red')
legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=pchs,col=colors,bty="n")

matplot(secTNF_sim[,1],secTNF_sim[,2:4],type='l',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (sim)',xlim=c(0,120),,ylim=c(0,2000),
        main = 'No sec regulation',lwd=2)#,yaxt='n')
#axis(2, at = c(0,1.2,3), labels=c(0,400,1000))

matplot(secTNF_sim_sec[,1],secTNF_sim_sec[,2:4],type='l',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (sim)',xlim=c(0,120),ylim=c(0,2000),
        main = 'Sec regulation only',lwd=2)#,yaxt ='n')
#axis(2, at = c(0,1.2,3), labels=c(0,400,1000))

matplot(secTNF_exp[,1],secTNF_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (exp)',xlim=c(0,120),lwd=2)
title(main = "Experimental data",col.main='red')
legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=pchs,col=colors,bty="n")

# tl only



#matplot(secTNF_sim_tl[,1],secTNF_sim_tl[,2:4],type='l',pch=pchs,col=colors,
#        lty=pchs,xlab='Time (mins)',ylab='secTNF (simulation)',xlim=c(0,120),ylim=c(0,3.5),
#        main = 'Tl regulation only')

# sec only
#matplot(proTNF_sim_sec[,1],proTNF_sim_sec[,2:4],type='l',pch=pchs,col=colors,
#        lty=pchs,xlab='Time (mins)',ylab='proTNF (simulation)',xlim=c(0,120),ylim=c(0,100),
#        main = 'Sec regulation only')


# both

#matplot(proTNF_sim_both[,1],proTNF_sim_both[,2:4],type='l',pch=pchs,col=colors,
#        lty=pchs,xlab='Time (mins)',ylab='proTNF (simulation)',xlim=c(0,120),ylim=c(0,100),
#        main = 'Tl and sec regulation')

#matplot(secTNF_sim_both[,1],secTNF_sim_both[,2:4],type='l',pch=pchs,col=colors,
#        lty=pchs,xlab='Time (mins)',ylab='secTNF (simulation)',xlim=c(0,120),ylim=c(0,3.5),
#        main = 'Tl and sec regulation')

# exp


# ending 


dev.off()
system('open fig4.pdf')

