function residues = calScoreCustom23(input_pars)
addpath('../../src/')
% expdata
nfkb_exp = csvread('../../expdata/nfkb.csv',1,0);

% params
plot_flag = 1;


pars = getParams(); % wt parameters
input_pars = 10.^input_pars;
pars('V_tr') = 1;%input_pars(6); 
pars('Km_tr') = input_pars(1);
pars('k_pr') = input_pars(2);
pr_fold_mko = input_pars(3);
pr_fold_tko = input_pars(4);
pars('n') = input_pars(5); 

kdeg = [.02 .02 .07]; % wt, mko, tko 
k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold_mko pars('k_pr')/pr_fold_tko]; 

% init     
yinit_all = zeros(2,3);
yinit_all(2,:) = yinit_all(1,:).*k_pr_all/kdeg(3);    
yinit_all(1,:) = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                  2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;
% time 
times = 0:.1:120;%nascent_all(:,1);

%% simulations
% wt
[t,wt]= ode15s(@ode23,times,yinit_all(:,1),[],[],nfkb_exp(:,1:2), ...
                       pars);

% mko 
pars('k_pr') = k_pr_all(2);
pars('Km_tr') = pars('Km_tr') *2; 
pars('kdeg_m') = kdeg(2);
[~,mko]= ode15s(@ode23,times,yinit_all(:,2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('Km_tr') = input_pars(1); 
pars('kdeg_m') = kdeg(3);
[~,tko]= ode15s(@ode23,times,yinit_all(:,3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);

% get date
simData_nascent = [wt(:,1) mko(:,1) tko(:,1)];
simData_mRNA = [wt(:,2) mko(:,2) tko(:,2)];

if plot_flag
    subplot(2,2,1)
    plot(times,simData_nascent)
    hold off
    subplot(2,2,2)
    plot(times,simData_mRNA)
    hold off
end


%% calculate score 
% residue1, peak time 
residues = zeros(1,16); % 8 features. 

% features of nascent profile
[max_val, max_ind] = max(simData_nascent);
max_time = (max_ind-1)*0.1;

residues(1) = (max_time(1) - 30)/5 ; % wt peak time 
residues(2) = (max_time(2) - 60)/10 ; % mko peak time 
residues(3) = (max_time(3) - 30)/5 ; % tko peak time 

residues(4) = (max_val(1) / max_val(2) - 0.68)/0.11; % peak_tko / peak_wt 
residues(5) = (max_val(1) / max_val(3) - 0.75)/0.22; % peak_mko /peak_tko 

%residues(6) = (simData_nascent(601,1)/mean([simData_nascent(601,2),simData_nascent(601,3)])- ...
%               .83)/.18;
residues(6) = (simData_nascent(601,1)/simData_nascent(601,3)-1.55)/ ...
    .67;
residues(7) = (simData_nascent(601,1)/simData_nascent(601,2)-0.59)/.07;
residues(8) = (simData_nascent(1201,1)/simData_nascent(1201,2)-.52)/.21;
residues(9) = (simData_nascent(1201,1)/simData_nascent(1201,3)-1.49)/.43;

% features of mRNA profile
[max_val, max_ind] = max(simData_mRNA);
max_time = (max_ind-1)*0.1;

residues(10) = (max_time(1) - 60)/10 ; % wt peak time 
residues(11) = (max_time(2) - 60)/10 ; % mko peak time 
residues(12) = (max_time(3) - 60)/10 ; % tko peak time 

residues(13) = (max_val(3) / max_val(1) - 0.73)/0.05; % peak_tko / peak_wt 
residues(14) = (max_val(2) / max_val(3) - 0.5)/0.06; % peak_mko /peak_tko 
residues(15) = (simData_mRNA(1201,1)/simData_mRNA(1201,2) - 1.98)/.23;
residues(16) = (simData_mRNA(1201,1)/simData_mRNA(1201,3) - 1.17)/.14;

%
inds = find(abs(residues)>1);
residues(inds) = ones(size(inds))*999;


