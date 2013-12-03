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

pdf(file='fig2.pdf', height=11, width=8.5, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12) 

layout(matrix(c(0,1,0,2,3,4,0,5,6),3,3,byrow=T))

# plot 1: input nkfb
matplot(nfkb_exp[,1],nfkb_exp[,2:4],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='NFkB ',xlim=c(0,120),main='Experimental data')

# plot 2: wt simulation
matplot(nascent_predict[,1],nascent_predict[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120))

# plot 3: less wt induction
matplot(nascent_predict_less_wt_induction[,1],nascent_predict_less_wt_induction[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120),main='1.5 fold less wt induction')

# plot 4: higher wt process
matplot(nascent_predict_higher_wt_process[,1],nascent_predict_higher_wt_process[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120),main='1.5 fold higher wt process')

# plot 5: Output exp
matplot(nascent_exp[,1],nascent_exp[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent ',xlim=c(0,120),main='Experimental data')

# plot 6: less mko process
matplot(nascent_mko_1.5_fold_less_process[,1],nascent_mko_1.5_fold_less_process[,2:4],type='l',col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent (sim) ',xlim=c(0,120),main='1.5 fold less mko process')


dev.off()
system('open fig2.pdf')