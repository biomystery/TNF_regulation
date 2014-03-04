####
# read 
####

# read simulation data 
nfkb_lps_exp <- read.csv('../expdata/nfkb_LPS_w_wo_feedback.csv')
nfkb_cpg_exp <- read.csv('../expdata/nfkb_CpG_w_wo_feedback.csv')
nfkb_lps_exp_std <- read.csv('../expdata/nfkb_LPS_w_wo_feedback_std.csv')
nfkb_cpg_exp_std <- read.csv('../expdata/nfkb_CpG_w_wo_feedback_std.csv')

nfkb_sim_nofeedback <- read.csv('./simData/nfkb_sim_nofeedback.csv',header=F)
nfkb_sim_feedback <- read.csv('./simData/nfkb_sim_feedback.csv',header=F)

####  
# plot 
####

pdf(file='Fig6.pdf', height=4, width=4, onefile=TRUE, family='Helvetica', paper='letter', pointsize=8) 

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))#, respect = TRUE)
ylim = c(0,1)
xlim = c(0,8)
xat = seq(0,8,2)
colors = c('steelblue4','darkolivegreen3','coral3')
pchs <- c(18,17,15)

##
#simulation plot
## 
# LPS plot
matplot(nfkb_sim_feedback[,1]/60,cbind(nfkb_sim_feedback[,2],nfkb_sim_nofeedback[,2])/max(nfkb_sim_feedback[,2]),
        type='l',pch=rep(pchs[1],2),col=rep(colors[1],2),lwd=2,
        lty=c(1,2),xlab='Time (hr)',ylab='NFkBn (sim)',xlim=xlim,ylim=ylim,xaxt ='n')
axis(1,at=xat,labels=xat)
title(main = "LPS stimulation (sim)")

# CpG plot
matplot(nfkb_sim_feedback[,1]/60,cbind(nfkb_sim_feedback[,3],nfkb_sim_nofeedback[,3])/max(nfkb_sim_feedback[,3]),
        type='l',pch=rep(1,3),col=rep(colors[2],2),lwd=2,
        lty=c(1,2),xlab='Time (hr)',ylab='NFkBn (sim)',xlim=xlim,ylim=ylim,xaxt ='n')
axis(1,at=xat)
title(main = "CpG stimulation (sim)")



##
#exp plot
## 
# LPS plot

matplot(nfkb_lps_exp[,1],nfkb_lps_exp[,2:3]/max(nfkb_lps_exp[,2:3]),
        type='o',pch=rep(pchs[1],2),col=rep(colors[1],2),lwd=2,
        lty=c(1,2),xlab='Time (hr)',ylab='NFkBn (exp)',xlim=xlim,ylim=ylim,xaxt ='n')
axis(1,at=xat)
ub = (nfkb_lps_exp_std + nfkb_lps_exp)/max(nfkb_lps_exp[,2:3])
lb = -nfkb_lps_exp_std + nfkb_lps_exp/max(nfkb_lps_exp[,2:3])
arrows(nfkb_lps_exp_std[,1],ub[,2], nfkb_lps_exp_std[,1], lb[,2], angle=90, code=3, length=.02,ylim=ylim)
arrows(nfkb_lps_exp_std[,1],ub[,3], nfkb_lps_exp_std[,1], lb[,3], angle=90, code=3, length=.02,ylim=ylim)

title(main = "LPS stimulation (exp)",col.main='red')
legend("bottom",c("wt","tnf ko"),lty=c(1,2),pch=rep(pchs[1],2),col=rep(colors[1],2),bty="n")

# CpG plot

matplot(nfkb_cpg_exp[,1],nfkb_cpg_exp[,2:3]/max(nfkb_cpg_exp[,2:3]),
        type='o',pch=rep(pchs[1],2),col=rep(colors[2],2),lwd=2,
        lty=c(1,2),xlab='Time (hr)',ylab='NFkBn (exp) ',xlim=xlim,ylim=ylim,xaxt ='n')
axis(1,at=xat)
ub = (nfkb_cpg_exp_std + nfkb_cpg_exp)/max(nfkb_cpg_exp[,2:3])
lb = -nfkb_cpg_exp_std + nfkb_cpg_exp/max(nfkb_cpg_exp[,2:3])
arrows(nfkb_cpg_exp_std[,1],ub[,2], nfkb_cpg_exp_std[,1], lb[,2], angle=90, code=3, length=.02,ylim=ylim)
arrows(nfkb_cpg_exp_std[,1],ub[,3], nfkb_cpg_exp_std[,1], lb[,3], angle=90, code=3, length=.02,ylim=ylim)

title(main = "CpG stimulation (exp)",col.main='red')
legend("bottom",c("wt","tnf ko"),lty=c(1,2),pch=rep(pchs[1],2),col=rep(colors[2],2),bty="n")

# PIC plot
#matplot(nfkb_sim_feedback[,1],cbind(nfkb_sim_feedback[,4],nfkb_sim_nofeedback[,4]),
#        type='l',pch=rep(1,3),col=rep(colors[3],2),lwd=2,
#        lty=c(1,2),xlab='Time (hr)',ylab='NFkBn (sim)',xlim=xlim,xaxt ='n')
#axis(1,at=xat)
#title(main = "PIC stimulation (sim)")

# end and open the pdf 
dev.off()
system('open Fig6.pdf')
