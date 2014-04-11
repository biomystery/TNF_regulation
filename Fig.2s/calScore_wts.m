function residues = calScore_wts(input_pars,nfkb_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters

pars('k_pr') = input_pars(1); % input parameters
pars('V_tr') = input_pars(2); % input parameters

yinit = pars('V_tr')* nfkb_exp(1,[2 4])/pars('k_pr');

times = 0:.1:120;%nascent_all(:,1);

%% simulations

% wt
[t,wt]= ode15s(@ode2s,times,yinit(1),[],[],nfkb_exp(:,1:2), ...
                       pars);


[t,tko]= ode15s(@ode2s,times,yinit(2),[],[],nfkb_exp(:,[1 4]), ...
                       pars);

% get date
simData = [wt tko];
simData = simData/max(simData(:,1)); % normalized 

if plot_flag
    figure
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
residues = abs(simData - expData(:,[2,4]));
residues = residues(:);

