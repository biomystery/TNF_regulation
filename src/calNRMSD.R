CalNRMSD <- function(simData, expData, sigmaTh)
    # calculate Normalized RMSD between simulation and exp. time courses data
    # http://en.wikipedia.org/wiki/Root-mean-square_deviation#Normalized_root-mean-square_deviation
    #
    # NRMSD = RMSD/(2*sigma), sigma >= sigmaTh (percent of mean)
    #
    # input:
    #      simData --- rows: wt, mko, tko; cols: timepoint1, timepoint2...
    #      expData --- rows: wt, wt_std, mko, mko_std, tko, tko_std; cols same as simData

    
    
