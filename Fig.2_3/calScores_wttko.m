function residues = calScores_wttko(input_pars,nfkb_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters
pars('Km_tr') =input_pars(1);% input_pars(2);%input_pars(2); % input parameters
pars('V_tr') = 1;%input_pars(1); 
pr_fold_mko =4;
km_fold_mko = 1;%input_pars(1);%input_pars(1);

pr_fold_tko = input_pars(2);
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


pars('k_pr') = k_pr_all(3);
pars('Km_tr') = k_m_all(3);
[t,tko]= ode15s(@ode2,times,yinit(3),[],[],nfkb_exp(:,[1 4]), ...
                       pars);

% get date
simData = [wt tko];
simData = simData/max(simData(:,1)); % normalized 

if plot_flag
    set(gcf,'DefaultAxesColorOrder',[54 100 139; 79 148 205]/255)    
    plot(times,simData,'linewidth',1.5)


    hold on ; 
    errorbar(expData(:,1),expData(:,2),expData(:,3),'o','color',[54 ...
                        100 139]/255)
    errorbar(expData(:,1),expData(:,6),expData(:,7),'*','color',[79 148 205]/255)    
    hold off
end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,[2 6]));
residues = residues(:);
residues(end-2) =  residues(end-2)*2;

