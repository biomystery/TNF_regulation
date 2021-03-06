% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
expData = csvread('../expdata/mRNA.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 
nascent_exp(:,2:end)= nascent_exp(:,2:end)/100;
% plot options 
plot_flag = 1; 

% pars
pars = getParams(); % wt parameters
                    %km_s = [pars('Km_tr') pars('Km_tr') pars('Km_tr')];
km_s = [pars('Km_tr') pars('Km_tr') pars('Km_tr')];
k_pr_all = [pars('k_pr') pars('k_pr')/4.5 pars('k_pr')/1.5]; 
kdeg_m = [.02 .02 .07]; % wt, mko, tko 

% initial conditions
yinit = nascent_exp(1,2:end) .* k_pr_all/.07;
times = 0:.1:120;%nascent_all(:,1);

% wt 
pars('k_pr') = k_pr_all(1);
pars('kdeg_m') = kdeg_m(1);
pars('kdeg_m_basal') = 0.07;

[t,wt]= ode15s(@ode3,times,yinit(:,1),[],[],nascent_exp(:,1:2), ...
               pars);


% mko 
pars('k_pr') = k_pr_all(2);
pars('kdeg_m') = kdeg_m(2);

[~,mko]= ode15s(@ode3,times,yinit(:,2),[],[],nascent_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('kdeg_m') = kdeg_m(3);

[~,tko]= ode15s(@ode3,times,yinit(:,3),[],[],nascent_exp(:,[1 ...
                    4]),pars);

% write the data into the file 
csvwrite('./simData/mRNA_all.csv',[t wt mko tko])

figure
plot(t,[wt(:,1) mko(:,1) tko(:,1)]/wt(1))
xlim([0 120])


%%
% run R
!R CMD BATCH fig3.R
!rm *.Ro* 
