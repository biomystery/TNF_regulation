nascent_predict_different_pr <- read.csv('./simData/different_pr.csv',header=F,
                                              col.names=c('Time_mins','wt','mko','tko'))
nascent_predict<-read.csv('./simData/same_pr.csv',header=F,
                          col.names=c('Time_mins','wt','mko','tko'))
nascent_exp <- read.csv('../expdata/nascent.csv')
nfkb_exp <- read.csv('../expdata/nfkb.csv')

nrmsd <- read.csv('./simData/nrmsd.csv',head = F)


############################################################
# plot fig2s.pdf
library(fields)

pdf(file='fig2s.pdf', height=14, width=8, onefile=TRUE, family='Helvetica', paper='letter', pointsize=18)

prFold <- seq(0.1,3,length.out = length(nrmsd))
image.plot(x=prFold,y=prFold,z=as.matrix(nrmsd),xlab='Fold reduction (mko)',ylab='Fold reduction (tko)',
           legend.lab = 'nRMSD', horizontal = TRUE)
dev.off()
system('open fig2s.pdf')
