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
mRNA_predict_wt <- read.csv('mRNA_all_wt.csv',header=F,
                         col.names=c('Time_mins','wt','mko','tko'))
mRNA_predict_tko <- read.csv('mRNA_all_tko.csv',header=F,
                         col.names=c('Time_mins','wt','mko','tko'))


# draw input data 
#par(mfrow=c(2,3))
pdf(file='fig3_wt_mko.pdf', height=6, width=12, onefile=TRUE,  paper='letter', pointsize=14) 

layout(matrix(c(1,2,3,0,4,0),2,3,byrow=T))
colors <- c('steelblue4','skyblue1','steelblue3')
pchs <- c(18,17,15)
yats <- c(0,25,50,75,100)
# Output:wt, 1
ymax= max(cbind(mRNA_predict_wt[,2:4],mRNA_predict[,2:4],mRNA_predict_tko[,2:4]))
ymax = 100

matplot(mRNA_predict_wt[,1],mRNA_predict_wt[,c(2,3,4)],type="l",xlab='Time (mins)'
     ,ylab='mRNA (sim)',main='Assume all mRNAs are stable',ylim=c(0,ymax),
        lty=rep(1,3),col= colors, lwd =2, yaxt ='n' )

axis(2,at = yats, label = yats*2)


# Output: tko, 2 


matplot(mRNA_predict_tko[,1],mRNA_predict_tko[,c(2,3,4)],type="l",xlab='Time (mins)'
        ,ylab='mRNA (sim)',main='Assume all mRNAs are unstable',ylim=c(0,ymax),
        lty=rep(1,3),col= colors, lwd =2 , yaxt ='n')
axis(2,at = yats, label = yats*2)
# Output: mRNA, 3 
matplot(mRNA_predict[,1],mRNA_predict[,c(2,3,4)],type="l",xlab='Time (mins)'
        ,ylab='mRNA (sim)',main='Assume stabilzation by Trif',ylim=c(0,ymax),
        lty=rep(1,3),col= colors, lwd =2, yaxt ='n' )
axis(2,at = yats, label = yats*2)


#plot()

# mRNA exp data
mRNA_expdata <- read.csv('../expdata/mRNA.csv')

matplot(mRNA_expdata[,1],mRNA_expdata[,c(2,4,6)],type="o",pch=pchs,col=colors,xlab='Time (mins)'
     ,lty=rep(1,3),ylab='mRNA (exp)',main='Experimental data',col.main = 'red',
        xlim=c(0,120),ylim=c(0,200),lwd=2)

legend("topleft",c("wt","mko","tko"),lty=rep(1,3),pch=pchs,col=colors,bty="n")

library(plotrix)
#ref : https://stat.ethz.ch/pipermail/r-help/2009-April/195287.html
plotCI(x=rep(mRNA_expdata[,1],3),y=as.vector(mRNA_expdata[,c(2,4,6)]),
       uiw=as.vector(mRNA_expdata[,c(3,5,7)]),
       col=rep(colors,each = nrow(mRNA_expdata)),add=T)
#dispersion(rep(mRNA_expdata[,1],3),as.vector(mRNA_expdata[,c(2,4,6)]),as.vector(mRNA_expdata[,c(3,5,7)]))
# nascnet
#matplot(nascent[,1],nascent[,c(2,4,6)],type='o',pch=pchs,col=colors,
#        lty=rep(1,3),xlab='Time (mins)',ylab='Nascent ',xlim=c(0,120),lwd=2,
#        main='Experimental data (input)')
# save as pdf
#dev.copy2pdf(file="fig3_wt_mko.pdf", width = 10, height = 7)


dev.off()
system('open fig3_wt_mko.pdf')

#end
#par(mfrow=c(1,1))

