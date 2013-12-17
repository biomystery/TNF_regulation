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

pdf(file='Fig6.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=14) 

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))#, respect = TRUE)

xlim = c(0,480)
xat = seq(0,480,60)
colors = c('steelblue4','darkolivegreen3','coral3')
pchs <- c(18,17,15)

##
#simulation plot
## 
# LPS plot
matplot(nfkb_sim_feedback[,1],cbind(nfkb_sim_feedback[,2],nfkb_sim_nofeedback[,2]),
        type='l',pch=rep(pchs[1],2),col=rep(colors[1],2),lwd=2,
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn (sim)',xlim=xlim,xaxt ='n')
axis(1,at=xat,labels=xat)
title(main = "LPS stimulation (sim)")

# CpG plot
matplot(nfkb_sim_feedback[,1],cbind(nfkb_sim_feedback[,3],nfkb_sim_nofeedback[,3]),
        type='l',pch=rep(1,3),col=rep(colors[2],2),lwd=2,
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn (sim)',xlim=xlim,xaxt ='n')
axis(1,at=xat)
title(main = "CpG stimulation (sim)")



##
#exp plot
## 
# LPS plot
nfkb_lps_exp_fold = t(rbind(nfkb_lps_exp[,2]/nfkb_lps_exp[1,2],nfkb_lps_exp[,3]/nfkb_lps_exp[1,3]))
matplot(nfkb_lps_exp[,1],nfkb_lps_exp_fold,
        type='o',pch=rep(pchs[1],2),col=rep(colors[1],2),lwd=2,
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn (exp)',xlim=xlim,xaxt ='n')
axis(1,at=xat)
title(main = "LPS stimulation (exp)",col.main='red')
legend("bottom",c("wt","tnf ko"),lty=c(1,2),pch=rep(pchs[1],2),col=rep(colors[1],2),bty="n")

# CpG plot
nfkb_cpg_exp_fold = t(rbind(nfkb_cpg_exp[,2]/nfkb_cpg_exp[1,2],nfkb_cpg_exp[,3]/nfkb_cpg_exp[1,3]))

matplot(nfkb_cpg_exp[,1],nfkb_cpg_exp_fold,
        type='o',pch=rep(pchs[1],2),col=rep(colors[2],2),lwd=2,
        lty=c(1,2),xlab='Time (mins)',ylab='NFkBn (exp) ',xlim=xlim,xaxt ='n')
axis(1,at=xat)
title(main = "CpG stimulation (exp)",col.main='red')
legend("bottom",c("wt","tnf ko"),lty=c(1,2),pch=rep(pchs[1],2),col=rep(colors[2],2),bty="n")


# end and open the pdf 
dev.off()
system('open Fig6.pdf')
