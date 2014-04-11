####
# read 
####

# read experimental data 
nfkb_exp <- read.csv('../expdata/nfkb.csv')
nascent_exp <- read.csv('../expdata/nascent.csv')
mRNA_exp <- read.csv('../expdata/mRNA.csv')
sec_exp <- read.csv('../expdata/TNF_secrection.csv')


# read simulation data
nfkb_sim_nf<- read.csv('./simData/nfkb_sim_nf.csv')
nascent_sim_nf <- read.csv('./simData/nascent_sim_nf.csv')
mRNA_sim_nf <- read.csv('./simData/mRNA_sim_nf.csv')
sec_sim_nf <- read.csv('./simData/sec_sim_nf.csv')

nfkb_sim <- read.csv('./simData/nfkb_sim.csv')
nascent_sim <- read.csv('./simData/nascent_sim.csv')
mRNA_sim <- read.csv('./simData/mRNA_sim.csv')
sec_sim <- read.csv('./simData/sec_sim.csv')


####  
# plot 
####

pdf(file='Fig5b.pdf', height=12, width=8, onefile=TRUE, family='Helvetica', paper='letter', pointsize=14) 

layout(matrix(c(5,1,6,2,7,3,8,4), 4, 2, byrow = TRUE))#, FFFR)

xlim = c(0,120)
colors <- c('steelblue4','skyblue1','steelblue3')
pchs <- c(18,17,15)
##
#exp plot
## 
# nfkb plot
xlim2 = c(0,240)
matplot(nfkb_exp[,1],nfkb_exp[,2:4],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB (exp)',xlim=xlim,lwd=2)

title(main = "Experimental data",col.main='red')

# nascnet plot
matplot(nascent_exp[,1],nascent_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (exp)',xlim=xlim,lwd=2)

# mRNA plot
matplot(mRNA_exp[,1],mRNA_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (exp)',xlim=xlim,lwd=2)


# sec exp plot 
matplot(sec_exp[,1],sec_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (exp)',xlim=xlim,lwd=2)
legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=pchs,col=colors,bty="n")

##
#sim plot
## 
# nfkb plot
matplot(nfkb_sim[,1],nfkb_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB (sim)',xlim=xlim2,lwd=2)
matlines(nfkb_sim_nf[,1],nfkb_sim_nf[,2:4],type='l',col=colors,
        lty=rep(2,3),xlab='Time (mins)',ylab='NFkB (sim)',xlim=xlim2,lwd=2)

title(main = "Simulation",col.main='red')

# nascnet plot
matplot(nascent_sim[,1],nascent_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim)',xlim=xlim2,lwd=2)
matlines(nascent_sim_nf[,1],nascent_sim_nf[,2:4],type='l',col=colors,
        lty=rep(2,3),xlab='Time (mins)',ylab='Nascent (sim)',xlim=xlim2,lwd=2)

# mRNA plot
matplot(mRNA_sim[,1],mRNA_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='mRNA (sim)',xlim=xlim2,lwd=2)
matlines(mRNA_sim_nf[,1],mRNA_sim_nf[,2:4],type='l',col=colors,
        lty=rep(2,3),xlab='Time (mins)',ylab='mRNA (sim)',xlim=xlim2,lwd=2)

#axis(2,at = c(0,.7,1.4)*.8,labels = c(0,40,80))
# sec exp plot 
matplot(sec_sim[,1],sec_sim[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='secTNF (sim)',xlim=xlim2,lwd=2)
legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=colors,bty="n")

matlines(sec_sim_nf[,1],sec_sim_nf[,2:4],type='l',col=colors,
        lty=rep(2,3),xlab='Time (mins)',ylab='secTNF (sim)',xlim=xlim2,lwd=2)

#axis(2,at = c(0,16,32,48),labels = c(0,400,800,1200))

# end and open the pdf 
dev.off()
system('open Fig5b.pdf')
