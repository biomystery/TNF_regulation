function residues = calScoreCustom(input_pars)

% expdata
nfkb_exp = csvread('../../expdata/nfkb.csv',1,0);
expData =  csvread('../../expdata/nascent2.csv',1,0);

% params
plot_flag = 1;
pars = getParams(); % wt parameters

input_pars = 10.^input_pars;

pars('V_tr') = input_pars(1);
pars('Km_tr') = input_pars(2);
pars('k_pr') = input_pars(3);
pr_fold_mko = input_pars(4);
pr_fold_tko = input_pars(5);
scale = input_pars(6)

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
[~,nascent_mko]= ode15s(@ode2,times,yinit_all(2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
[~,nascent_tko]= ode15s(@ode2,times,yinit_all(3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);

% get date
simData = [nascent_wt nascent_mko nascent_tko];

if plot_flag
    subplot(2,2,1)
    %plot(expData(:,1),expData(:,[2,4,6]),'*')
    %hold on 
    plot(times,simData)
    hold off
end

%simData = simData((expData(:,1))*10+1,:);
%simData = simData;

%% calculate score 
% residue1, peak time 
residues = zeros(1,9); % 8 features. 
simData = simData;
[max_val, max_ind] = max(simData);
max_time = (max_ind-1)*0.1;

residues(1) = (max_time(1) - 30)/10 ; % wt peak time 
residues(2) = (max_time(2) - 60)/20 ; % mko peak time 
residues(3) = (max_time(3) - 30)/10 ; % tko peak time 

residues(4) = (max_val(3) / max_val(1) - 1.25)/0.1; % peak_wt / peak_mko 
residues(5) = (max_val(2) / max_val(3) - 1.05)/0.01; % peak_mko /peak_tko 
residues(6) = (simData(601,1)/mean([simData(601,2),simData(601,3)])- ...
               .75)/.15;
residues(7) = (simData(1201,1)/mean([simData(1201,2),simData(1201,3)])-.75)/.15;
residues(8) = (simData(1201,1)-simData(601,1) + 6)/1;
residues(9) = (max_val(3) -25)/2; 