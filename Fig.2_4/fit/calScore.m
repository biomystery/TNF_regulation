function rmsd_residues = calScore(input_pars)
addpath('../../src/')
% expdata
nfkb_exp = csvread('../../expdata/nfkb.csv',1,0);
nascent_exp = csvread('../../expdata/nascent.csv',1,0);
mRNA_exp = csvread('../../expdata/mRNA.csv',1,0);
secTNF_exp = csvread('../../expdata/TNF_secrection.csv',1,0);
proTNF_exp = csvread('../../expdata/combineTNF.csv',1,0);


mRNA_exp = mRNA_exp(1:4,:);

mRNA_exp(:,2:end) = mRNA_exp(:,2:end)/max(mRNA_exp(:,2));
nascent_exp(:,2:end) = nascent_exp(:,2:end)/max(nascent_exp(:,2));
secTNF_exp(:,2:end) = secTNF_exp(:,2:end)/max(secTNF_exp(:,2));
secTNF_exp = secTNF_exp(1:4,:);
proTNF_exp(:,2:end) = proTNF_exp(:,2:end)/max(proTNF_exp(:,2));
proTNF_exp = [proTNF_exp(:,1:2) proTNF_exp(:,2)*.2 proTNF_exp(:,3) proTNF_exp(:,3)*.2 ...
           proTNF_exp(:,4) proTNF_exp(:,4)*.2]; 

% params
plot_flag = 1;

pars = getParams(); % wt parameters

input_pars = 10.^input_pars;
pars('V_tr') = 1;%input_pars(6); 
pars('Km_tr') = input_pars(1);
pars('k_pr') = input_pars(2);
pars('pr_fold_mko') = input_pars(3);
pars('pr_fold_tko') = input_pars(4);
pars('n') = 2; 
pars('Km_tr_fold') = input_pars(5);

pars('k_tl') = input_pars(6);
pars('tl_fold') = input_pars(7);
pars('sec_fold') = input_pars(8);
pars('k_sec') = input_pars(9);
pars('kdeg_P') = input_pars(10);

k_pr_all = [pars('k_pr') pars('k_pr')/(pars('pr_fold_mko')) pars('k_pr')/pars('pr_fold_tko')]; 
kdeg = [.02 .02 .07]; % wt, mko, tko 
k_tls = [pars('k_tl'),pars('k_tl'),pars('k_tl')/pars('tl_fold')]; 

yinit_all = zeros(2,3);
yinit_all(1,:) = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                  2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;
yinit_all(2,:) = yinit_all(1,:).*k_pr_all/kdeg(3);
yinit_all(3,:) = yinit_all(2,:).*k_tls/pars('kdeg_p');

% time 
times = 0:.1:120;%nascent_all(:,1);

%% simulations
% wt 
[t,wt]= ode15s(@ode2_4,times,yinit_all(:,1),[],[],nfkb_exp(:,1:2), ...
               pars);
% mko 
pars('k_pr') = k_pr_all(2);
back_km = pars('Km_tr');
pars('Km_tr') = pars('Km_tr') *pars('Km_tr_fold'); 
pars('kdeg_m') = kdeg(2);
[~,mko]= ode15s(@ode2_4,times,yinit_all(:,2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('kdeg_m') = kdeg(3);
pars('Km_tr') = back_km; 
pars('k_tl') = k_tls(3); % set as tko
[~,tko]= ode15s(@ode2_4,times,yinit_all(:,3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);


% get date
simData_nascent = [wt(:,1) mko(:,1) tko(:,1)]/max(wt(:,1));
simData_mRNA = [wt(:,2) mko(:,2) tko(:,2)]/max(wt(:,2));
simData_proTNF = [wt(:,end) mko(:,end) tko(:,end)]/max(wt(:,end));

%% secretion 
pars = getParams(); % set as default
k_secs = [pars('k_sec'),pars('k_sec'),pars('k_sec')/ ...
          pars('sec_fold')];

yinit_all(3,:) = yinit_all(2,:).*k_tls./(pars('kdeg_p')+k_secs); % init 

% wt 
pars('k_pr') = k_pr_all(1);
pars('kdeg_m') = kdeg(1);
pars('Km_tr') = back_km; 
pars('k_tl') = k_tls(1);

[t,wt]= ode15s(@ode2_4a,times,yinit_all(:,1),[],[],nfkb_exp(:,1:2), ...
               pars);
% mko 
pars('k_pr') = k_pr_all(2);
back_km = pars('Km_tr');
pars('Km_tr') = pars('Km_tr') *pars('Km_tr_fold'); 
pars('kdeg_m') = kdeg(2);
[~,mko]= ode15s(@ode2_4a,times,yinit_all(:,2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('kdeg_m') = kdeg(3);
pars('Km_tr') = back_km; 

pars('k_sec') = k_secs(3); % set as tko
pars('k_tl') = k_tls(3); % set as tko
[~,tko]= ode15s(@ode2_4a,times,yinit_all(:,3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);


simData_secTNF = cumsum([wt(:,end)*k_secs(1) mko(:,end)*k_secs(2) ...
                    tko(:,end)*k_secs(3)]);
simData_secTNF = simData_secTNF/max(simData_secTNF(:,1));



% plot
if plot_flag
    subplot(2,2,1)
    errorbar(nascent_exp(:,1),nascent_exp(:,2),nascent_exp(:,3),'o','Color',[54 100 139]/255)
    hold on 
    errorbar(nascent_exp(:,1),nascent_exp(:,4),nascent_exp(:,5),'^','Color',[135 206 255]/255)
    errorbar(nascent_exp(:,1),nascent_exp(:,6),nascent_exp(:,7),'*','Color',[79 148 205]/255)

    plot(times,simData_nascent)
    hold off;
    xlim([0 120])
    ylim([0 2.0])
    xlabel('Time (min)'); ylabel('Nascent mRNA (a.u.)')
    
    subplot(2,2,2)
    errorbar(mRNA_exp(:,1),mRNA_exp(:,2),mRNA_exp(:,3),'o','Color',[54 100 139]/255)
    hold on 
    errorbar(mRNA_exp(:,1),mRNA_exp(:,4),mRNA_exp(:,5),'^','Color',[135 206 255]/255)
    errorbar(mRNA_exp(:,1),mRNA_exp(:,6),mRNA_exp(:,7),'*','Color',[79 148 205]/255)

    plot(times,simData_mRNA)
    hold off
    xlim([0 120])
    xlabel('Time (min)'); ylabel('Mature mRNA (a.u.)')
    subplot(2,2,3)
    errorbar(proTNF_exp(:,1),proTNF_exp(:,2),proTNF_exp(:,3),'o','Color',[54 100 139]/255)
    hold on 
    errorbar(proTNF_exp(:,1),proTNF_exp(:,4),proTNF_exp(:,5),'^','Color',[135 206 255]/255)
    errorbar(proTNF_exp(:,1),proTNF_exp(:,6),proTNF_exp(:,7),'*','Color',[79 148 205]/255)

    plot(times,simData_proTNF)
    hold off
    xlim([0 120])
    xlabel('Time (min)'); ylabel('ProTNF (a.u.)')
    subplot(2,2,4)
    errorbar(secTNF_exp(:,1),secTNF_exp(:,2),secTNF_exp(:,3),'o','Color',[54 100 139]/255)
    hold on 
    errorbar(secTNF_exp(:,1),secTNF_exp(:,4),secTNF_exp(:,5),'^','Color',[135 206 255]/255)
    errorbar(secTNF_exp(:,1),secTNF_exp(:,6),secTNF_exp(:,7),'*','Color',[79 148 205]/255)

    plot(times,simData_secTNF)
    xlabel('Time (min)'); ylabel('secTNF (a.u.)')
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

simData_proTNF = simData_proTNF((proTNF_exp(:,1))*10+1,:);
[rmsd_residues_proTNF, nrmsd_residues_proTNF] = calNRMSDResidues(simData_proTNF,proTNF_exp(:, ...
                                                  2:end));


simData_secTNF = simData_secTNF((secTNF_exp(:,1))*10+1,:);
[rmsd_residues_secTNF, nrmsd_residues_secTNF] = calNRMSDResidues(simData_secTNF,secTNF_exp(:, ...
                                                  2:end));

rmsd_residues = [rmsd_residues_nascent; rmsd_residues_mRNA;rmsd_residues_proTNF;rmsd_residues_secTNF];
