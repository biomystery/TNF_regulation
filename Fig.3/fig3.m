% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
expData = csvread('../expdata/mRNA.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 

% plot options 
plot_flag = 1; 

% pars
pars = getParams(); % wt parameters
pars('k_pr') = 0.9;
k_pr_all = [pars('k_pr') pars('k_pr')/1 pars('k_pr')/1]; 
kdeg_m = [.02 .02 .07]; % wt, mko, tko 
k_pr = pars('k_pr');
kdeg_m_mko = [0.07, 0.02];

% initial conditions
yinit = nascent_exp(1,2:end) .* k_pr_all/kdeg_m(3);
times = 0:.1:120;%nascent_all(:,1);

% wt 
[t,wt]= ode15s(@ode3,times,yinit(:,1),[],[],nascent_exp(:,1:2), ...
               pars);
pars('k_pr') = k_pr_all(1);
pars('kdeg_m') = kdeg_m(1);

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
%csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])

figure
plot(t,[wt(:,1) mko(:,1) tko(:,1)]/wt(1))
hold on 
plot(expData(:,1),expData(:,[2,4,6]),'o')
xlim([0 120])
hold off
