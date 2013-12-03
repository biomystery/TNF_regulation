##
# To show roles of transcriptional regulation and
# RNA stabilization
##

# parameters
interp_method = 'linear'

# input data 
nascent <- read.csv('nascent.csv')
halflife <- read.csv('halflife.csv')#,row.names=1)
kdeg_all <- log(2)/halflife

mRNA_predict <- read.csv('mRNA_all2.csv',header=F,
                         col.names=c('Time_mins','wt','mko','tko'))
mRNA_predict_wt <- read.csv('mRNA_all_wt.csv',header=F,
                         col.names=c('Time_mins','wt','mko','tko'))
mRNA_predict_tko <- read.csv('mRNA_all_tko.csv',header=F,
                         col.names=c('Time_mins','wt','mko','tko'))


# draw input data 
#par(mfrow=c(2,3))
pdf(file='fig3_wt_mko.pdf', height=10, width=8.5, onefile=TRUE, family='Helvetica', paper='letter', pointsize=14) 

layout(matrix(c(0,5,0,1,2,3,0,4,0),3,3,byrow=T))

# Output:wt, 1
ymax= max(cbind(mRNA_predict_wt[,2:4],mRNA_predict[,2:4],mRNA_predict_tko[,2:4]))

plot(mRNA_predict_wt$Time_mins,mRNA_predict_wt$wt,type="l",xlab='Time (mins)'
     ,ylab='mRNA',main='Assume all mRNAs are stable',ylim=c(0,ymax))
lines(mRNA_predict_wt$Time_mins,mRNA_predict_wt$mko,nascent$mko,type="l",col='purple')
lines(mRNA_predict_wt$Time_mins,mRNA_predict_wt$tko,type="l",col='cyan3')

mtext('mRNA\'= nascent(t) - kdeg* mRNA'
      ,adj=0,cex=.55,col='grey')
legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=rep(NA,3),col=c('black','purple','cyan3'),bty="n")

# Output: tko, 2 
plot(mRNA_predict_tko$Time_mins,mRNA_predict_tko$wt,type="l",xlab='Time (mins)'
     ,ylab='mRNA',main='Assume all mRNAs are unstable',ylim=c(0,ymax))
lines(mRNA_predict_tko$Time_mins,mRNA_predict_tko$mko,nascent$mko,type="l",col='purple')
lines(mRNA_predict_tko$Time_mins,mRNA_predict_tko$tko,type="l",col='cyan3')

mtext('mRNA\'= nascent(t) - kdeg* mRNA'
      ,adj=0,cex=.55,col='grey')

# Output: mRNA, 3 
plot(mRNA_predict$Time_mins,mRNA_predict$wt,type="l",xlab='Time (mins)'
     ,ylab='mRNA',main='Assume stabilzation by Trif',ylim=c(0,ymax))
lines(mRNA_predict$Time_mins,mRNA_predict$mko,nascent$mko,type="l",col='purple')
lines(mRNA_predict$Time_mins,mRNA_predict$tko,type="l",col='cyan3')

mtext('mRNA\'= nascent(t) - kdeg* mRNA'
      ,adj=0,cex=.55,col='grey')


#plot()

# mRNA exp data
mRNA_expdata <- read.csv('mRNA.csv')

plot(mRNA_expdata$Time_mins,mRNA_expdata$wt_fold,type="b",xlab='Time (mins)'
     ,ylab='mRNA',main='Experimental data (output)',xlim=c(0,120),ylim=c(0,200))
lines(mRNA_expdata$Time_mins,mRNA_expdata$mko_fold,nascent$mko,type="b",col='purple')
lines(mRNA_expdata$Time_mins,mRNA_expdata$tko_fold,type="b",col='cyan3')

mtext('black:wt, purple:mko, cyan:tko'
      ,adj=0,cex=.55,col='grey')


# nascnet
matplot(nascent[,1],nascent[,c(2,4,6)],type='b',pch=rep(1,3),col=c('black','purple','cyan3'),
        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent ',xlim=c(0,120),main='Experimental data (input)')
# save as pdf
#dev.copy2pdf(file="fig3_wt_mko.pdf", width = 10, height = 7)


dev.off()
system('open fig3_wt_mko.pdf')

#end
#par(mfrow=c(1,1))

