% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData =  csvread('../expdata/nascent.csv',1,0);

plot_flag = 0; 
pars = getParams(); % wt parameters
pr_fold = [pars('pr_fold') 1];
k_pr = pars('k_pr');



pr_fold = linspace(0.1,3,N);
rmsd = zeros(N,N);
nrmsd = zeros(N,N);


k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold(i) pars('k_pr')/pr_fold(j)]; 
yinit_all = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                  2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;

times = 0:.1:120;%nascent_all(:,1);

% wt
[t,nascent_wt]= ode15s(@ode2,times,yinit_all(1),[],[],nfkb_exp(:,1:2), ...
                       pars);

% mko 
pars('k_pr') = k_pr_all(2);
[~,nascent_mko]= ode15s(@ode2,times,yinit_all(2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
[~,nascent_tko]= ode15s(@ode2,times,yinit_all(3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);

pars('k_pr') = k_pr;

% get date
simData = [nascent_wt nascent_mko nascent_tko];
simData = simData((expData(:,1))*10+1,:);

% calculate score 
[rmsd, nrmsd] = calNRMSD(simData,expData(:,2:end));


