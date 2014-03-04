function rmsd_residues = calScore(input_pars)

% expdata
nfkb_exp = csvread('../../expdata/nfkb.csv',1,0);
expData =  csvread('../../expdata/nascent.csv',1,0);

% params
plot_flag = 1;
pars = getParams(); % wt parameters

input_pars = 10.^input_pars;

pars('V_tr') = 1;%input_pars(1);

pars('Km_tr') = input_pars(1);
pars('k_pr') = input_pars(2);
pr_fold_mko = input_pars(3);
pr_fold_tko = input_pars(4);
pars('n') = input_pars(5); 
scale = input_pars(6);
%scale = input_pars(6)

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

% rescaling:
simData = simData/scale;


% plot 
if plot_flag
    %    subplot(2,2,1)
    errorbar(expData(:,1),expData(:,2),expData(:,3),'*','Color',[0 0 1])
    hold on 
    errorbar(expData(:,1),expData(:,4),expData(:,5),'*','Color',[0 .5 ...
                        0])
    errorbar(expData(:,1),expData(:,6),expData(:,7),'*','Color',[1 0 ...
                        0 ])
    plot(times,simData)
    hold off
    xlim([0 120])
    ylim([0 30])
end


% calculate score 
simData = simData((expData(:,1))*10+1,:);
[rmsd_residues, nrmsd_residues] = calNRMSDResidues(simData,expData(:, ...
                                                  2:end));




