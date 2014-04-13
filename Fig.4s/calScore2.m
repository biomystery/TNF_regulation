function residues = calScore2(input_pars,mRNA_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters

pars('kdeg_p') = input_pars(1); % input parameters
pars('k_tl') = 1; 



fold_tl_tko =input_pars(2); % input parameters
fold_sec_tko = 1;%input_pars(2);

ktls = [pars('k_tl') pars('k_tl') pars('k_tl')/fold_tl_tko];


yinit = mRNA_exp(1,2:end).*ktls./(pars('kdeg_p'));
times = 0:.1:120;%nascent_all(:,1);

%% simulations
% wt
[t,wt]= ode15s(@ode4,times,yinit(:,1),[],[],mRNA_exp(:,1:2), ...
                       pars);

[t,mko]= ode15s(@ode4,times,yinit(:,2),[],[],mRNA_exp(:,[1,3]), ...
                       pars);
pars('k_tl') = ktls(3); 

[t,tko]= ode15s(@ode4,times,yinit(:,3),[],[],mRNA_exp(:,[1,4]), ...
                       pars);

% get date

simData = [wt mko tko];
simData = simData/max(wt); % normalized 

if plot_flag
    plot(times,simData,'linewidth',1.5)
    hold on ; 
    plot(expData(:,1),expData(:,2),'o','color',[54 100 139]/255)
    plot(expData(:,1),expData(:,3),'^','color',[135 206 255]/255)
    plot(expData(:,1),expData(:,4),'*','color',[79 148 205]/255)
    hold off; 
end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,2:end));
residues = residues(:); 

