%% figs1. The figure for same and different processing rate 

% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData = csvread('../expdata/nascent.csv',1,0);
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

plot_flag = 1; 
pars = getParams(); % wt parameters
pr_fold = [1.5 4.5]; 
k_pr = pars('k_pr');

N = 10;
M = 19;
pr_rate = [linspace(0.01,.1,10) linspace(0.2,1,9)];
km_rate = linspace(0.1,1,N);

score = zeros(M,N);

for i =1:M % different pr vs same pr. 
    for j = 1:N
    input_pars(1) = pr_rate(i);
    input_pars(2) = km_rate(j);
    residues = calScore_wt(input_pars,nascent_exp,expData,0);
    score(i,j) = sum(residues.^2);
    end
    disp([i,j])
end



save ./simData/fig2s_wt.mat 
csvwrite('./simData/fig2s_wt.csv',score)

% write the data into the file 
%csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])
figure
if plot_flag
    imagesc(km_rate,log10(pr_rate),score)
    xlabel('k_{M}');ylabel('log10(k_{Pr})');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 
i = find(score == min(score)) ; 

calScore_wt(pr_rate(i),nascent_exp,expData,1);


%% run R 

% run R
%!R CMD BATCH fig3s2.r
%!rm *.Ro* 

