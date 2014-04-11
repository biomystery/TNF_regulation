nascent_predict <- read.csv('./simData/fig2s2.csv',header=F,
                                              col.names=c('Time_mins','wt','mko','tko'))
nascent_predictb <- read.csv('./simData/fig2s2b.csv',header=F,
                                              col.names=c('Time_mins','wt','mko','tko'))

nascent_exp <- read.csv('../expdata/nascent.csv')
nfkb_exp <- read.csv('../expdata/nfkb.csv')

nrmsd <- read.csv('./simData/nrmsd.csv')


############################################################
# plot fig2.pdf

pdf(file='fig2s2.pdf', height=3, width=6, onefile=F, family='Helvetica', paper='special', pointsize=10) 
layout(matrix(c(1,2),1,2,byrow=T))

colors <- c('steelblue4','skyblue1','steelblue3')
pchs <- c(18,17,15)
ylims <- c(0,ceiling(max(nascent_exp[,-1]))+5)

# plot 2: wt simulation
matplot(nascent_predict[,1],nascent_predict[,2:4],type='l',col=colors,lwd=2,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120),
        main='fold_mko = 4.5, fold_tko =1.5')

matplot(nascent_predictb[,1],nascent_predictb[,2:4],type='l',col=colors,lwd=2,
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120),
        main='with km_tr 2* in mko')

legend('bottom',c("wt","mko","tko"),lty=rep(1,3),col=colors,bty="n",cex=.8)


dev.off()

system('open fig2s2.pdf')

