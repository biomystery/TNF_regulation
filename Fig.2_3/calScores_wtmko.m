function residues = calScores_wtmko(input_pars,nfkb_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters
pars('V_tr') = 0.7; 
pars('Km_tr') =.5;% input_pars(2);%input_pars(2); % input parameters

pr_fold_mko =input_pars(1);
km_fold_mko =1;%input_pars(2);%input_pars(1);%input_pars(1);

pr_fold_tko = 1;
km_fold_tko= 1;%input_pars(2);%input_pars(2);

k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold_mko pars('k_pr')/pr_fold_tko]; 
k_m_all = [pars('Km_tr') pars('Km_tr')*km_fold_mko pars('Km_tr')/ ...
           km_fold_tko];

yinit = pars('V_tr')* nfkb_exp(1,2:4).^pars('n')./(nfkb_exp(1,2:4).^pars('n')+k_m_all.^pars('n'))./k_pr_all;


times = 0:.1:120;%nascent_all(:,1);

%% simulations

% wt
pars('k_pr') = k_pr_all(1);
[t,wt]= ode15s(@ode2,times,yinit(1),[],[],nfkb_exp(:,1:2), ...
                       pars);

pars('k_pr') = k_pr_all(2);
pars('Km_tr') = k_m_all(2);
[t,mko]= ode15s(@ode2,times,yinit(2),[],[],nfkb_exp(:,[1 3]), ...
                       pars);


% get date
simData = [wt mko];
simData = simData/max(simData(:,1)); % normalized 

if plot_flag
    plot(times,simData,'linewidth',1.5)

    hold on ; 
% $$$     errorbar(expData(:,1),expData(:,2),expData(:,3),'o','color',[54 ...
% $$$                         100 139]/255)
% $$$     errorbar(expData(:,1),expData(:,4),expData(:,5),'^','color',[135 206 255]/255)
% $$$     hold off
end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,[2 4]));
residues = residues(:);


