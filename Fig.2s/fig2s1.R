nascent_predict_different_pr <- read.csv('./simData/different_pr.csv',header=F,
                                              col.names=c('Time_mins','wt','mko','tko'))
nascent_predict<-read.csv('./simData/same_pr.csv',header=F,
                          col.names=c('Time_mins','wt','mko','tko'))
nascent_exp <- read.csv('../expdata/nascent.csv')
nfkb_exp <- read.csv('../expdata/nfkb.csv')

nrmsd <- read.csv('./simData/nrmsd.csv')


############################################################
# plot fig2.pdf

pdf(file='fig2s1.pdf', height=8, width=6, onefile=F, family='Helvetica', paper='special', pointsize=10) 

layout(matrix(c(0,1,1,0,2,2,3,3,0,4,4,0),3,4,byrow=T))
colors <- c('steelblue4','skyblue1','steelblue3')
pchs <- c(18,17,15)
ylims <- c(0,ceiling(max(nascent_exp[,-1]))+5)
# plot 1: input nkfb
matplot(nfkb_exp[,1],nfkb_exp[,2:4],type='o',pch=pchs,col=colors,lwd=2,
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB ',xlim=c(0,120),
        main='Experimental data (input)')

# plot 2: wt simulation
matplot(nascent_predict[,1],nascent_predict[,2:4],type='l',col=colors,lwd=2,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120),
        main='Same RNA process in knockouts & wt')

#legend("bottomright",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=c('black','purple','cyan3'),bty="n")

# plot 3: lower process in knockouts. 
matplot(nascent_predict_different_pr[,1],nascent_predict_different_pr[,2:4],type='l',col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',
        xlim=c(0,120),lwd=2,main='1.5 fold less process in knockouts')
# legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=colors,bty="n")

# plot 4: Output exp
matplot(nascent_exp[,1],nascent_exp[,c(2,4,6)],type='o',pch=pchs,col=colors,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (Exp)',xlim=c(0,120),
        ylim= ylims,lwd=2)#,main='Experimental data (Output)')
legend('bottom',c("wt","mko","tko"),lty=rep(1,3),pch=pchs,col=colors,bty="n",cex=.8)


dev.off()

system('open fig2s1.pdf') 
