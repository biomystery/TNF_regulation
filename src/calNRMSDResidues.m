function  [rmsd_residues,nrmsd_residues] = calNRMSDResidues(simData,expData,sigmaRatioTh)
% calculate Normalized RMSD between simulation and exp. time courses data
% http://en.wikipedia.org/wiki/Root-mean-square_deviation%Normalized_root-mean-square_deviation
%
% NRMSD = RMSD/(2*sigma), sigma >= sigmaTh (percent of mean)
%
% input:
%      simData --- rows: wt, mko, tko; cols: timepoint1, timepoint2...
%      expData --- rows: wt, wt_std, mko, mko_std, tko, tko_std;
%      cols same as simData

    % default value for sigmaRatioTh
    if nargin <3
        sigmaRatioTh = 0.05;
    end

    expData_mean = expData(:,[1,3,5]);
    expData_sigma = expData(:,[2,4,6]);
    expData_sigmaRatio = expData_sigma./expData_mean;

    % adjusted by the sigma threshold
    expData_sigma(expData_sigmaRatio<=sigmaRatioTh) = ...
        expData_mean(expData_sigmaRatio<=sigmaRatioTh) * sigmaRatioTh;


    rmsd_residues = simData(:) - expData_mean(:);
    nrmsd_residues = (simData - expData_mean)./expData_sigma;
    nrmsd_residues = nrmsd_residues(:);
end