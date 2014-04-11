function residues = calScoreb(input_pars,mRNA_exp,expData,plot_flag)

% params
%plot_flag = 1;
pars = getParams(); % wt parameters

pars('kdeg_p') = input_pars(1); % input parameters
pars('k_tl') = 1; 
pars('k_sec') = input_pars(2);
pars('km_sec') = 1; 

fold_tl_tko =1;%input_pars(2); % input parameters
fold_sec_tko =1;% input_pars(2);

ktls = [pars('k_tl') pars('k_tl') pars('k_tl')/fold_tl_tko];
ksecs = [pars('k_sec') pars('k_sec') pars('k_sec')/fold_sec_tko];

yinit = mRNA_exp(1,2:end).*ktls./(pars('kdeg_p') + ksecs);%pars('k_sec')); % init 
times = 0:.1:120;%nascent_all(:,1);

%% simulations
% wt
[t,wt]= ode15s(@ode4_1,times,yinit(:,1),[],[],mRNA_exp(:,1:2), ...
                       pars);

[t,mko]= ode15s(@ode4_1,times,yinit(:,2),[],[],mRNA_exp(:,[1,3]), ...
                       pars);
pars('k_tl') = ktls(3); 
pars('k_secs') = ksecs(3); 
[t,tko]= ode15s(@ode4_1,times,yinit(:,3),[],[],mRNA_exp(:,[1,4]), ...
                       pars);

% get date
wt = cumsum(wt./(wt+pars('km_sec')))*ksecs(1);%pars('k_sec');
mko = cumsum(mko./(mko+pars('km_sec')))*ksecs(2);%pars('k_sec');
tko = cumsum(tko./(tko+pars('km_sec')))*ksecs(3);%pars('k_sec');

simData = [wt mko tko];
simData = simData/max(wt); % normalized 

if plot_flag
    plot(times,simData,'linewidth',1.5)
    hold on ; 
    errorbar(expData(:,1),expData(:,2),expData(:,3),'o','color',[54 100 139]/255)
    errorbar(expData(:,1),expData(:,4),expData(:,5),'^','color',[135 206 255]/255)
    errorbar(expData(:,1),expData(:,6),expData(:,7),'*','color',[79 148 205]/255)
    hold off; 
end

%% calculate score 
% residue1, peak time 
simData = simData(expData(:,1)*10+1,:);
residues = abs(simData - expData(:,[2,4,6]));
residues = residues(:); 

