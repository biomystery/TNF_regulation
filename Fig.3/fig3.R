##
# To show roles of transcriptional regulation and
# RNA stabilization
##

# parameters
interp_method = 'linear'

# input data
nascent2 <- read.csv('nascent2.csv')
halflife <- read.csv('../expdata/halflife.csv')#,row.names=1)
kdeg_all <- log(2)/halflife

mRNA_predict <- read.csv('mRNA_all.csv',header=F,
                         col.names=c('Time_mins','wt','mko','tko'))

# draw input data 
par(mfrow=c(2,2))

# input: nascent activity (w/wo interpretion,data from fig2.R)
plot(nascent2$x,nascent2$wt,type="l",col='black',xlab='Time (mins)'
     ,ylab='nascent',main='Input')
lines(nascent2$x,nascent2$mko,type="l",col='purple')
lines(nascent2$x,nascent2$tko,type="l",col='cyan3')
mtext('Predicted by interped NFkB data'
      ,adj=0,cex=.55,col='grey')

# Output: mRNA 
plot(mRNA_predict$Time_mins,mRNA_predict$wt,type="l",xlab='Time (mins)'
     ,ylab='mRNA',main='Onput')
lines(mRNA_predict$Time_mins,mRNA_predict$mko,type="l",col='purple')
lines(mRNA_predict$Time_mins,mRNA_predict$tko,type="l",col='cyan3')

#lines(nascent2$x,nascent2$wt,type="l",col='black')
#lines(nascent2$x,nascent2$mko,type="l",col='purple')
#lines(nascent2$x,nascent2$tko,type="l",col='cyan3')
mtext('mRNA\'= nascent(t) - kdeg* mRNA'
      ,adj=0,cex=.55,col='grey')

# mRNA exp data
mRNA_expdata <- read.csv('../expdata/mRNA.csv')
mRNA_expdata

plot(mRNA_expdata$Time_mins,mRNA_expdata$wt_fold,type="b",xlab='Time (mins)'
     ,ylab='mRNA',main='Experimental data',xlim=c(0,120))
lines(mRNA_expdata$Time_mins,mRNA_expdata$mko_fold,nascent$mko,type="b",col='purple')
lines(mRNA_expdata$Time_mins,mRNA_expdata$tko_fold,type="b",col='cyan3')
mtext('black:wt, purple:mko, cyan:tko'
      ,adj=0,cex=.55,col='grey')
# save as pdf
dev.copy2pdf(file="fig3.pdf", width = 8, height = 7)
system('open fig3.pdf')

#end
par(mfrow=c(1,1))

