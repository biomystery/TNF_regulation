function residues = calScoreCustom(input_pars,nascent_exp,expData)

% params
plot_flag = 0;
pars = getParams(); % wt parameters

pr_fold_mko = input_pars(1);
pr_fold_tko = input_pars(2);
kdeg_m = [.02 .02 .07]; % wt, mko, tko 


k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold_mko pars('k_pr')/pr_fold_tko]; 
yinit = nascent_exp(1,2:end) .* k_pr_all/kdeg_m(3);

times = 0:.1:120;%nascent_all(:,1);

%% simulations

% wt
pars('k_pr') = k_pr_all(1);
pars('kdeg_m') = kdeg_m(1);
[t,wt]= ode15s(@ode3,times,yinit(:,1),[],[],nascent_exp(:,1:2), ...
                       pars);

% mko 
pars('k_pr') = k_pr_all(2);
pars('kdeg_m') = kdeg_m(2);
[~,mko]= ode15s(@ode3,times,yinit(:,2),[],[],nascent_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('kdeg_m') = kdeg_m(3);
[~,tko]= ode15s(@ode3,times,yinit(:,3),[],[],nascent_exp(:,[1 ...
                    4]),pars);

% get date
simData = [wt mko tko];

if plot_flag
    figure
    subplot(2,2,1)
    %plot(expData(:,1),expData(:,[2,4,6]),'*')
    %hold on 
    plot(times,simData)
    hold off
end

%% calculate score 
% residue1, peak time 
residues = zeros(1,7); % 7 features. 
[max_val, max_ind] = max(simData);
max_time = (max_ind-1)*0.1;

residues(1) = (max_time(1) - 60)/10 ; % wt peak time 
residues(2) = (max_time(2) - 60)/10 ; % mko peak time 
residues(3) = (max_time(3) - 60)/10 ; % tko peak time 
residues(4) = (max_val(3) / max_val(1) - 0.73)/0.05; % peak_tko / peak_wt 
residues(5) = (max_val(2) / max_val(3) - 0.5)/0.06; % peak_mko /peak_tko 
residues(6) = (simData(1201,1)/simData(1201,2) - 1.98)/.23;
residues(7) = (simData(1201,1)/simData(1201,3) - 1.17)/.14;


