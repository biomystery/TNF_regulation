####
# read 
####

# read experimental data 
nfkb_exp <- read.csv('../expdata/nfkb.csv')
nascent_exp <- read.csv('../expdata/nascent.csv')
mRNA_exp <- read.csv('../expdata/mRNA.csv')
sec_exp <- read.csv('../expdata/TNF_secrection.csv')


# read simulation data 
nfkb_sim <- read.csv('./simData/nfkb_sim.csv')
nascent_sim <- read.csv('./simData/nascent_sim.csv')
mRNA_sim <- read.csv('./simData/mRNA_sim.csv')
sec_sim <- read.csv('./simData/sec_sim.csv')


####  
# plot 
####

pdf(file='Fig5b_TLR4_LPS.pdf', height=8, width=8, onefile=TRUE, family='Helvetica', paper='letter', pointsize=14) 

layout(matrix(c(4,1,5,2,6,3), 3, 2, byrow = TRUE))#, FFFR)

xlim = c(0,120)
colors <- c('steelblue4','skyblue1','steelblue3')
pchs <- c(18,17,15)
##
#exp plot
## 
# nfkb plot
matplot(nfkb_exp[,1],nfkb_exp[,2:4],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB (exp)',xlim=xlim,lwd=2)
title(main = "Experimental data",col.main='red')

# nascnet plot
matplot(nascent_exp[,1],nascent_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (exp)',xlim=xlim,lwd=2)

# mRNA plot
#matplot(mRNA_exp[,1],mRNA_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
#        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (exp)',xlim=xlim,lwd=2)


# sec exp plot 
matplot(sec_exp[,1],sec_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (exp)',xlim=xlim,lwd=2)
legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=pchs,col=colors,bty="n")

##
#sim plot
## 
# nfkb plot
matplot(nfkb_sim[,1],nfkb_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB (sim)',xlim=xlim,lwd=2)

title(main = "Simulation",col.main='red')

# nascnet plot
matplot(nascent_sim[,1],nascent_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim)',xlim=xlim,lwd=2,yaxt = 'n')
axis(2,at = c(0,.2,.4,.6,.8,1.0),labels = c(0,5,10,15,20,25))

# mRNA plot
#matplot(mRNA_sim[,1],mRNA_sim[,2:4],type='l',col=colors,
#        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (sim)',xlim=xlim,lwd=2,yaxt = 'n')

#axis(2,at = c(0,.7,1.4)*.8,labels = c(0,40,80))
# sec exp plot 
matplot(sec_sim[,1],sec_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (sim)',xlim=xlim,lwd=2,yaxt = 'n')
legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=colors,bty="n")
axis(2,at = c(0,16,32,48),labels = c(0,400,800,1200))

# end and open the pdf 
dev.off()
system('open Fig5b_TLR4_LPS.pdf')
