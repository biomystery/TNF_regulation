% read data from R 
addpath('../src/')
simDataOrg = csvread('./simData/same_pr.csv',0,0);
expData =  csvread('../expdata/nascent.csv',1,0);
simData = simDataOrg((expData(:,1))*10+1,2:end);

% calculate score 
[rmsd, nrmsd] = calNRMSD(simData,expData(:,2:end));


% plot 
if plot_flag
    imagesc(pr_fold,pr_fold,nrmsd)
end 

% run R
%!R CMD BATCH fig2.R
%!rm *.Ro* 
%!rm *.Rh*
%!rm *.RD*