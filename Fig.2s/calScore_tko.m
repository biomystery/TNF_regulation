function residues = calScore_tko(input_pars,nfkb_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters
pars('n') = 2; 
pars('k_pr') = 0.4;
pars('Km_tr') = input_pars(2); % input parameters
k_pr_all = [pars('k_pr') pars('k_pr')/input_pars(1)];

yinit = pars('V_tr')* nfkb_exp(1,[2 4]).^pars('n')./(nfkb_exp(1,[2 4]).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;

times = 0:.1:120;%nascent_all(:,1);

%% simulations

% wt
[t,wt]= ode15s(@ode2,times,yinit(:,1),[],[],nfkb_exp(:,1:2), ...
                       pars);

pars('k_pr') = k_pr_all(2);

[t,tko]= ode15s(@ode2,times,yinit(:,2),[],[],nfkb_exp(:,[1 4]), ...
                       pars);

% get date
simData = [wt tko];
simData = simData/max(wt); % normalized 

if plot_flag
    set(gcf,'DefaultAxesColorOrder',[54 100 139; 79 148 205]/255)        
    plot(times,simData,'linewidth',1.5)

    hold on ; 

    errorbar(expData(:,1),expData(:,2),expData(:,3),'o','color',[54 ...
                        100 139]/255)
    
    errorbar(expData(:,1),expData(:,6),expData(:,7),'*','color',[79 ...
                        148 205]/255)

end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,[2 6])); % wt tko
residues = residues(:); 

