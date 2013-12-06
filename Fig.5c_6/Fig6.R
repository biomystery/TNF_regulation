####
# read 
####

# read simulation data 
nfkb_lps_exp <- read.csv('../expdata/nfkb_LPS_w_wo_feedback.csv')
nfkb_cpg_exp <- read.csv('../expdata/nfkb_CpG_w_wo_feedback.csv')

nfkb_sim_nofeedback <- read.csv('./simData/nfkb_sim_nofeedback.csv',header=F)
nfkb_sim_feedback <- read.csv('./simData/nfkb_sim_feedback.csv',header=F)

####  
# plot 
####

pdf(file='Fig6.pdf', height=6, width=12, onefile=TRUE, family='Helvetica', paper='letter', pointsize=14) 

layout(matrix(c(1,2,3,4,5,0), 2, 3, byrow = TRUE))#, respect = TRUE)

xlim = c(0,480)
xat = seq(0,480,60)

##
#simulation plot
## 
# LPS plot
matplot(nfkb_sim_feedback[,1],cbind(nfkb_sim_feedback[,2],nfkb_sim_nofeedback[,2]),
        type='l',pch=rep(1,3),col=c('black','red'),
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn ',xlim=xlim,xaxt ='n')
axis(1,at=xat)
title(main = "LPS stimulation")

# CpG plot
matplot(nfkb_sim_feedback[,1],cbind(nfkb_sim_feedback[,3],nfkb_sim_nofeedback[,3]),
        type='l',pch=rep(1,3),col=c('black','red'),
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn ',xlim=xlim,xaxt ='n')
axis(1,at=xat)
title(main = "CpG stimulation")

# PIC plot
matplot(nfkb_sim_feedback[,1],cbind(nfkb_sim_feedback[,4],nfkb_sim_nofeedback[,4]),
        type='l',pch=rep(1,3),col=c('black','red'),
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn ',xlim=xlim,xaxt ='n')
axis(1,at=xat)
title(main = "PIC stimulation")

legend("topright",c("w TNF feedback","wo TNF feedback"),lty=c(1,2),pch=rep(NA,2),col=c('black','red'),bty="n")

##
#exp plot
## 
# LPS plot
nfkb_lps_exp_fold = t(rbind(nfkb_lps_exp[,2]/nfkb_lps_exp[1,2],nfkb_lps_exp[,3]/nfkb_lps_exp[1,3]))
matplot(nfkb_lps_exp[,1],nfkb_lps_exp_fold,
        type='b',pch=rep(1,3),col=c('black','red'),
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn (fold) ',xlim=xlim,xaxt ='n')
axis(1,at=xat)
title(main = "LPS stimulation")

# CpG plot
nfkb_cpg_exp_fold = t(rbind(nfkb_cpg_exp[,2]/nfkb_cpg_exp[1,2],nfkb_cpg_exp[,3]/nfkb_cpg_exp[1,3]))

matplot(nfkb_cpg_exp[,1],nfkb_cpg_exp_fold,
        type='b',pch=rep(1,3),col=c('black','red'),
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn (fold) ',xlim=xlim,xaxt ='n')
axis(1,at=xat)
title(main = "CpG stimulation")


# end and open the pdf 
dev.off()
system('open Fig6.pdf')
