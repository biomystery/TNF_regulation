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

pdf(file='Fig5b_TLR4_LPS.pdf', height=9, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=14) 

layout(matrix(c(5,1,6,2,7,3,8,4), 4, 2, byrow = TRUE))#, respect = TRUE)

xlim = c(0,120)
##
#exp plot
## 
# nfkb plot
matplot(nfkb_exp[,1],nfkb_exp[,2:4],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB ',xlim=xlim)
title(main = "Experimental data",col.main='red')

# nascnet plot
matplot(nascent_exp[,1],nascent_exp[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent ',xlim=xlim)

# mRNA plot
matplot(mRNA_exp[,1],mRNA_exp[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA ',xlim=xlim)


# sec exp plot 
matplot(sec_exp[,1],sec_exp[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF',xlim=xlim)

##
#sim plot
## 
# nfkb plot
matplot(nfkb_sim[,1],nfkb_sim[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB ',xlim=xlim)
title(main = "Simulation",col.main='red')

# nascnet plot
matplot(nascent_sim[,1],nascent_sim[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent ',xlim=xlim)

# mRNA plot
matplot(mRNA_sim[,1],mRNA_sim[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA ',xlim=xlim)


# sec exp plot 
matplot(sec_sim[,1],sec_sim[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF',xlim=xlim)

# end and open the pdf 
dev.off()
system('open Fig5b_TLR4_LPS.pdf')
