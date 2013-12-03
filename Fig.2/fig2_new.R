nascent_predict_less_wt_induction <- read.csv('nascent_wt_1.5_fold_less_induction.csv',header=F,
                                              col.names=c('Time_mins','wt','mko','tko'))

nascent_predict_higher_wt_process <- read.csv('nascent_wt_1.5_fold_higher_process.csv',header=F,
                                              col.names=c('Time_mins','wt','mko','tko'))

nascent_mko_1.5_fold_less_process <- read.csv('nascent_mko_1.5_fold_less_process.csv',header=F,
                                              col.names=c('Time_mins','wt','mko','tko'))

nascent_predict<-read.csv('nascent_wt.csv',header=F,
                          col.names=c('Time_mins','wt','mko','tko'))

nascent_exp <- read.csv('../Fig.2/nascent.csv')
nfkb_exp <- read.csv('nfkb_input.csv')

pdf(file='fig2.pdf', height=11, width=8.5, onefile=TRUE, family='Helvetica', paper='letter', pointsize=18) 

layout(matrix(c(0,1,1,0,2,2,3,3,0,4,4,0),3,4,byrow=T))

# plot 1: input nkfb
matplot(nfkb_exp[,1],nfkb_exp[,2:4],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB ',xlim=c(0,120),main='Experimental data (input)')

# plot 2: wt simulation
matplot(nascent_predict[,1],nascent_predict[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120),main='Same RNA process in knockouts & wt')

legend("bottom",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=c('black','purple','cyan3'),bty="n")

# plot 4: higher wt process
matplot(nascent_predict_higher_wt_process[,1],nascent_predict_higher_wt_process[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120),main='1.5 fold less process in knockouts')

# plot 5: Output exp
matplot(nascent_exp[,1],nascent_exp[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent ',xlim=c(0,120),main='Experimental data (Output)')



dev.off()
system('open fig2.pdf')