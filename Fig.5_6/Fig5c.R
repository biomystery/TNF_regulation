####
# read 
####

# read experimental data 
mRNA_exp <- read.csv('../expdata/mRNA_LPS_CpG_PIC.csv')
sec_exp <- read.csv('../expdata/elisa_LPS_CpG_PIC.csv')


# read simulation data 
mRNA_sim_no_feedback <- read.csv('./simData/mRNA_sim.csv')
sec_sim_no_feedback <- read.csv('./simData/sec_sim.csv')
mRNA_sim <- read.csv('./simData/mRNA_sim_feedback.csv')
sec_sim <- read.csv('./simData/sec_sim_feedback.csv')


####  
# plot 
####

pdf(file='Fig5c.pdf', height=6, width=12, onefile=TRUE, family='Helvetica', paper='letter', pointsize=14) 

layout(matrix(c(3,1,4,2), 2, 2, byrow = TRUE))#, respect = TRUE)

xlim = c(0,240)
xat = seq(0,240,60)
colors = c('steelblue4','darkolivegreen3','coral3')
pchs <- c(18,17,15)
##
#exp plot
## 
# mRNA plot
matplot(mRNA_exp[,1],mRNA_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (exp)',xlim=xlim,xaxt ='n',lwd=2)
axis(1,at=xat)
title(main = "Experimental data",col.main='red')

# elisa plot
matplot(sec_exp[,1],sec_exp[,2:4],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Secreted TNF (exp)',xlim=xlim,xaxt ='n',lwd=2)
axis(1,at=xat)

legend("topleft",c("LPS","CpG","PIC"),lty=rep(1,3),pch=pchs,col=colors,bty="n")


##
#Sim plot
## 
# mRNA plot
matplot(mRNA_sim[,1],mRNA_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (sim)',xlim=xlim,xaxt ='n',lwd=2)
matlines(mRNA_sim_no_feedback[,1],mRNA_sim_no_feedback[,2:4],type='l',col=colors,
        lty=rep(2,3),lwd=2)#,xlim=xlim,xaxt ='n',lwd=2)

axis(1,at=xat)
title(main = "Simulation")

# elisa plot
matplot(sec_sim[,1],sec_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Secreted TNF (sim)',xlim=xlim,xaxt ='n',lwd=2,ylim=c(0,1.5))
matlines(sec_sim_no_feedback[,1],sec_sim_no_feedback[,2:4],type='l',col=colors,
        lty=rep(2,3),lwd=2)#,xlab='Time (mins)',ylab='Secreted TNF (sim)',xlim=xlim,xaxt ='n',lwd=2)

axis(1,at=xat)
legend("topleft",c("LPS","CpG","PIC"),lty=rep(1,3),pch=rep(NA,3),col=colors,bty="n")


# end and open the pdf 
dev.off()
system('open Fig5c.pdf')
