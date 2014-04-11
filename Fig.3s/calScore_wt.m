function residues = calScore_wt(input_pars,nascent_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters
pars('k_pr') = input_pars(1); % input parameters
yinit = nascent_exp(1,2) * pars('k_pr')/0.07; % init

times = 0:.1:120;%nascent_all(:,1);

%% simulations

% wt
[t,wt]= ode15s(@ode3,times,yinit(:,1),[],[],nascent_exp(:,1:2), ...
                       pars);
% get date
simData = wt;
simData = simData/max(simData); % normalized 

if plot_flag
    %    figure
    errorbar(expData(:,1),expData(:,2),expData(:,3),'ob')
    hold on ; 
    plot(times,simData)
    hold off
    csvwrite('./simData/best_fit_wt.csv',[times' simData])
    csvwrite('./simData/exp_fit_wt.csv',expData)    
end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,2));

