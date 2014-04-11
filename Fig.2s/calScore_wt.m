function residues = calScore_wt(input_pars,nfkb_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters
pars('n') = 2; 
pars('k_pr') = input_pars(1); % input parameters
pars('Km_tr') = input_pars(2); % input parameters

yinit = pars('V_tr')* nfkb_exp(1,2)^pars('n')./(nfkb_exp(1,2)^pars('n')+pars('Km_tr')^pars('n'))/pars('k_pr');

times = 0:.1:120;%nascent_all(:,1);

%% simulations

% wt
[t,wt]= ode15s(@ode2,times,yinit(:,1),[],[],nfkb_exp(:,1:2), ...
                       pars);
% get date
simData = wt;
simData = simData/max(simData); % normalized 

if plot_flag
    plot(times,simData,'linewidth',1.5)
    hold on ; 

    errorbar(expData(:,1),expData(:,2),expData(:,3),'o','color',[54 ...
                        100 139]/255)
end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,2));

