function rmsd_residues= calScore(input_pars)
addpath('../../src/')
% expdata
nfkb_exp = csvread('../../expdata/nfkb.csv',1,0);
nascent_exp =  csvread('../../expdata/nascent.csv',1,0);
mRNA_exp =  csvread('../../expdata/mRNA.csv',1,0);
mRNA_exp = mRNA_exp(1:4,:);
mRNA_exp(:,2:end) = mRNA_exp(:,2:end)/max(max(mRNA_exp(:,2:end)));
nascent_exp(:,2:end) = nascent_exp(:,2:end)/max(max(nascent_exp(:,2:end)));

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
nascent_scale = input_pars(6);
mRNA_scale = input_pars(7);
pars('Km_tr_fold') = input_pars(8);

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
pars('Km_tr') = pars('Km_tr') * pars('Km_tr_fold'); 
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
simData_nascent = simData_nascent/nascent_scale;
simData_mRNA = [wt(:,2) mko(:,2) tko(:,2)];
simData_mRNA = simData_mRNA/mRNA_scale;

if plot_flag
    subplot(2,2,1)
    errorbar(nascent_exp(:,1),nascent_exp(:,2),nascent_exp(:,3),'*','Color',[0 0 1])
    hold on 
    errorbar(nascent_exp(:,1),nascent_exp(:,4),nascent_exp(:,5),'*','Color',[0 .5 ...
                        0])
    errorbar(nascent_exp(:,1),nascent_exp(:,6),nascent_exp(:,7),'*','Color',[1 0 ...
                        0 ])

    plot(times,simData_nascent)
    
    hold off
    xlim([0 120])
    ylim([0 1])
    %
    subplot(2,2,2)
    errorbar(mRNA_exp(:,1),mRNA_exp(:,2),mRNA_exp(:,3),'*','Color',[0 0 1])
    hold on 
    errorbar(mRNA_exp(:,1),mRNA_exp(:,4),mRNA_exp(:,5),'*','Color',[0 .5 ...
                        0])
    errorbar(mRNA_exp(:,1),mRNA_exp(:,6),mRNA_exp(:,7),'*','Color',[1 0 ...
                        0 ])

    plot(times,simData_mRNA)
    hold off
    xlim([0 120])
end


%% calculate score 
simData_nascent = simData_nascent((nascent_exp(:,1))*10+1,:);
[rmsd_residues_nascent, nrmsd_residues_nascent] = calNRMSDResidues(simData_nascent,nascent_exp(:, ...
                                                  2:end));

simData_mRNA = simData_mRNA((mRNA_exp(1:4,1))*10+1,:);
[rmsd_residues_mRNA, nrmsd_residues_mRNA] = calNRMSDResidues(simData_mRNA,mRNA_exp(1:4, ...
                                                  2:end));

rmsd_residues = [rmsd_residues_nascent; rmsd_residues_mRNA];


