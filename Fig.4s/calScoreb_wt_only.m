function residues = calScoreb_wt_only(input_pars,mRNA_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters

pars('kdeg_p') = 1;%input_pars(1); % input parameters
pars('k_tl') = 1; 
pars('k_sec') = 100;%input_pars(2);
pars('km_sec') = input_pars(2); 


fold_tl_tko = 1;%input_pars(2); % input parameters
ktls = [pars('k_tl') pars('k_tl') pars('k_tl')/fold_tl_tko];

yinit = mRNA_exp(1,2:end).*ktls/(pars('kdeg_p') + pars('k_sec')); % init 
times = 0:.1:120;%nascent_all(:,1);

%% simulations
% wt
[t,wt]= ode15s(@ode4_1,times,yinit(:,1),[],[],mRNA_exp(:,1:2), ...
                       pars);


% get date
wt = cumsum(wt)*pars('k_sec');

simData = wt;
simData = simData/max(wt); % normalized 

if plot_flag
    plot(times,simData,'linewidth',1.5)
    hold on ; 
    errorbar(expData(:,1),expData(:,2),expData(:,3),'o','color',[54 100 139]/255)

end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,2));
residues(:) = residues(:); 

