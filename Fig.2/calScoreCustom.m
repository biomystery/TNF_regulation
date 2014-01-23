function residues = calScoreCustom(input_pars,nfkb_exp,expData)

% expdata
%nfkb_exp = csvread('../../expdata/nfkb.csv',1,0);
%expData =  csvread('../../expdata/nascent2.csv',1,0);

% params
plot_flag = 0;
pars = getParams(); % wt parameters

pr_fold_mko = 4.5;%input_pars(1);
pr_fold_tko = 1.5;%input_pars(2)
pars('Km_tr') = input_pars(1); 
km_back = pars('Km_tr');
pars('k_pr') = input_pars(2); 

k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold_mko pars('k_pr')/pr_fold_tko]; 
yinit_all = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                  2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;

times = 0:.1:120;%nascent_all(:,1);

%% simulations
% wt
[t,nascent_wt]= ode15s(@ode2,times,yinit_all(1),[],[],nfkb_exp(:,1:2), ...
                       pars);
% mko 
pars('k_pr') = k_pr_all(2);
pars('Km_tr') = pars('Km_tr') * 2; 

[~,nascent_mko]= ode15s(@ode2,times,yinit_all(2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('Km_tr') = km_back; 

[~,nascent_tko]= ode15s(@ode2,times,yinit_all(3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);

% get date
simData_nascent = [nascent_wt nascent_mko nascent_tko];

if plot_flag
    figure
    subplot(2,2,1)
    %plot(expData(:,1),expData(:,[2,4,6]),'*')
    %hold on 
    plot(times,simData_nascent)
    hold off
end

%% calculate score 
% residue1, peak time 
residues = zeros(1,9); % 9 features. 
[max_val, max_ind] = max(simData_nascent);
max_time = (max_ind-1)*0.1;

residues(1) = (max_time(1) - 30)/5 ; % wt peak time 
residues(2) = (max_time(2) - 60)/10 ; % mko peak time 
residues(3) = (max_time(3) - 30)/5 ; % tko peak time 
residues(4) = (max_val(1) / max_val(2) - 0.68)/0.11; % peak_tko / peak_wt 
residues(5) = (max_val(1) / max_val(3) - 0.75)/0.22; % peak_mko /peak_tko 
residues(6) = (simData_nascent(601,1)/simData_nascent(601,3)-1.55)/ ...
    .67;
residues(7) = (simData_nascent(601,1)/simData_nascent(601,2)-0.59)/.07;
residues(8) = (simData_nascent(1201,1)/simData_nascent(1201,2)-.52)/.21;
residues(9) = (simData_nascent(1201,1)/simData_nascent(1201,3)-1.49)/.43;
