function residues = calScore(input_pars,nascent_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters

pr_fold_mko = input_pars(1);%input_pars(1);
pr_fold_tko = input_pars(2);

kdeg_m = [.02 .02 .07]; % wt, mko, tko 
                        %pars('k_pr') = input_pars(1);
kdeg_m_basal = [0.07 0.07 0.07]; 

k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold_mko pars('k_pr')/pr_fold_tko]; 
yinit = nascent_exp(1,2:end) .* k_pr_all./kdeg_m_basal;


times = 0:.1:120;%nascent_all(:,1);

%% simulations

% wt
pars('k_pr') = k_pr_all(1);
pars('kdeg_m') = kdeg_m(1);
pars('kdeg_m_basal') = kdeg_m_basal(1);

[t,wt]= ode15s(@ode3,times,yinit(:,1),[],[],nascent_exp(:,1:2), ...
                       pars);

% mko 
 pars('k_pr') = k_pr_all(2);
 pars('kdeg_m') = kdeg_m(2);
 pars('kdeg_m_basal') = kdeg_m_basal(2);
 [~,mko]= ode15s(@ode3,times,yinit(:,2),[],[],nascent_exp(:,[1 ...
                     3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('kdeg_m') = kdeg_m(3);
pars('kdeg_m_basal') = kdeg_m_basal(3);
[~,tko]= ode15s(@ode3,times,yinit(:,3),[],[],nascent_exp(:,[1 ...
                    4]),pars);

% get date
simData = [wt mko tko];
simData = simData/max(simData(:,1));  % normalized

if plot_flag
    plot(times,simData,'linewidth',1.5)

    hold on 
    errorbar(expData(:,1),expData(:,2),expData(:,3),'o','color',[54 100 139]/255)
    errorbar(expData(:,1),expData(:,4),expData(:,5),'^','color',[135 206 255]/255)
    errorbar(expData(:,1),expData(:,6),expData(:,7),'*','color',[79 148 205]/255)

    hold off

    %    csvwrite('./simData/exp_fit.csv',expData)    
    %    set(gca,'markersize',10)

end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,[2,4,6]));
residues = residues(:);


