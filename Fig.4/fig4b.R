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
pdf(file='fig4b.pdf', height=4, width=8, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12) 
colors <- c('steelblue4','skyblue1','steelblue3')
pchs <- c(18,17,15)

layout(matrix(1:6,2,3,byrow=T))

# no regulation
genotypes <- cbind('wt','mko','tko')

tl_reg <- cbind(max(proTNF_sim[,4])/max(proTNF_sim[,2]), max(proTNF_sim_tl[,4])/max(proTNF_sim_tl[,2]),max(proTNF_exp[,4])/max(proTNF_exp[,2]))

barx= barplot(tl_reg,names.arg =
    c('No tl regulation','tl regulation','experimental data'),main =
    'trif ko / wt',ylab=('Max proTNF ratio'),col = 'steelblue3') 
arrows(barx[3],.3371*(1+0.25), barx[3], .3371*(1-0.25), angle=90, code=3, length=.2)

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


dev.off()
system('open fig4b.pdf')

