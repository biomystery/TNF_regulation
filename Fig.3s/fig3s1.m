% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
expData = csvread('../expdata/mRNA.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 

% plot options 
plot_flag = 0; 

% pars
pars = getParams(); % wt parameters
                    %km_s = [pars('Km_tr') pars('Km_tr')
                    %pars('Km_tr')];
pr_fold = [1.5 3 4.5];
for i = 1:3
k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold(i) pars('k_pr')/pr_fold(i)]; 
kdeg_m = [.02 .02 .07]; % wt, mko, tko 
k_pr = pars('k_pr');


% initial conditions
yinit = nascent_exp(1,2:end) .* k_pr_all/kdeg_m(3);
times = 0:.1:120;%nascent_all(:,1);

% wt 
pars('k_pr') = k_pr_all(1);
pars('kdeg_m') = kdeg_m(1);
[t,wt{i}]= ode15s(@ode3s,times,yinit(:,1),[],[],nascent_exp(:,1:2), ...
               pars);


% mko 
pars('k_pr') = k_pr_all(2);
pars('kdeg_m') = kdeg_m(2);
[~,mko{i}]= ode15s(@ode3s,times,yinit(:,2),[],[],nascent_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('kdeg_m') = kdeg_m(3);
[~,tko{i}]= ode15s(@ode3s,times,yinit(:,3),[],[],nascent_exp(:,[1 ...
                    4]),pars);
end

% write the data into the file 
csvwrite('./simData/mRNA_all.csv',[t wt{1} mko{1} tko{1}  mko{2} ...
                   tko{2}  mko{3} tko{3}])


%%
% run R
!R CMD BATCH fig3s1.R
!rm *.Ro* 
