% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 

% plot options 
plot_flag = 1; 

% pars
pars = getParams(); % wt parameters
k_pr_all = [pars('k_pr') pars('k_pr')/6 pars('k_pr')]; 
kdeg_m = [.02 .02 .07]; % wt, mko, tko 
k_pr = pars('k_pr');
kdeg_m_mko = [0.07, 0.02]

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
plot(t,wt(:,1),'k')
hold on 
plot(t,tko(:,1),'c')
plot(t,mko(:,1),'Color', [0.5 0 0.5])
hold off
