##
# To show roles of transcriptional regulation and
# RNA stabilization
##

# parameters
interp_method = 'linear'

# input data 
nascent <- read.csv('../expdata/nascent.csv')
halflife <- read.csv('../expdata/halflife.csv')#,row.names=1)
kdeg_all <- log(2)/halflife

mRNA_predict <- read.csv('mRNA_all2.csv',header=F,
                         col.names=c('Time_mins','wt','mko','tko'))
head(mRNA_predict)

# draw input data 
par(mfrow=c(2,2))

# input: nascent activity (w/wo interpretion,data from fig2.R)
plot(nascent$Time_mins,nascent$mko_fold,type="b",col='purple',xlab='Time (mins)'
     ,ylab='nascent',main='Input')
lines(nascent$Time_mins,nascent$wt_fold,type="b")
lines(nascent$Time_mins,nascent$tko_fold,type="b",col='cyan3')


mtext('Use nascent data directly'
      ,adj=0,cex=.55,col='grey')

# Output: mRNA 
plot(mRNA_predict$Time_mins,mRNA_predict$wt,type="l",xlab='Time (mins)'
     ,ylab='mRNA',main='Onput')
lines(mRNA_predict$Time_mins,mRNA_predict$mko,nascent$mko,type="l",col='purple')
lines(mRNA_predict$Time_mins,mRNA_predict$tko,type="l",col='cyan3')

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
dev.copy2pdf(file="fig3_2.pdf", width = 8, height = 7)
system('open fig3_2.pdf')

#end
par(mfrow=c(1,1))

