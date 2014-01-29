function residues = calScoreCustom2_4(input_pars)
addpath('../../src/')
% expdata
nfkb_exp = csvread('../../expdata/nfkb.csv',1,0);
nascent_exp = csvread('../../expdata/nascent.csv',1,0);
mRNA_exp = csvread('../../expdata/mRNA.csv',1,0);
sec_exp = csvread('../../expdata/TNF_secrection.csv',1,0);
pro_exp = csvread('../../expdata/combineTNF.csv',1,0);

% params
plot_flag = 1;


pars = getParams(); % wt parameters
pars('kdeg_p') = log(2)/15;
pars('k_sec') =  log(2)/15;


input_pars = 10.^input_pars;
ps = [ 0.65    0.4    4.2    1.5    2];
%ps = [ 0.6491    0.4230    4.1937    1.5414    2.0581];
pars('V_tr') = 1;%input_pars(6); 
pars('Km_tr') = ps(1);
pars('k_pr') = ps(2);
pr_fold_mko = ps(3);
pr_fold_tko = ps(4);
pars('n') = ps(5); 
%pars('Km_tr_fold') = input_pars(6);
pars('k_tl') = input_pars(1);
pars('tl_fold') = input_pars(2);
pars('sec_fold') = input_pars(3);
 pars('k_sec') = input_pars(4);
 pars('k_degP') = input_pars(5);

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
pars('Km_tr') = pars('Km_tr') *2; 
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
simData_nascent = [wt(:,1) mko(:,1) tko(:,1)];
simData_mRNA = [wt(:,2) mko(:,2) tko(:,2)];
simData_proTNF = [wt(:,end) mko(:,end) tko(:,end)];

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
pars('Km_tr') = pars('Km_tr') *2; 
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


simData_secTNF = cumsum([wt(:,end)*k_secs(1) mko(:,end)*k_secs(2) tko(:,end)*k_secs(3)]);


% plot
if plot_flag
    subplot(2,2,1)
    plot(times,simData_proTNF)
    hold off
    subplot(2,2,2)
    plot(times,simData_secTNF)
    hold off
end


%% calculate score 
% residue1, peak time 
residues = zeros(1,30); % 8 features. 

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
residues(16) = (simData_mRNA(1201,1)/simData_mRNA(1201,3) - 1.17)/ ...
    .14;

% feature of proTNF
[max_val, max_ind] = max(simData_proTNF);
max_time = (max_ind-1)*0.1;

residues(17) = (max_time(1) - 60)/10 ; % wt peak time 
residues(18) = (max_time(2) - 60)/10 ; % mko peak time 
residues(19) = (max_time(3) - 60)/10 ; % tko peak time 
residues(20) = (max_val(3) / max_val(1) - 0.34)/0.06; % peak_tko / peak_wt 
residues(21) = (max_val(2) / max_val(3) - 0.82)/0.14; % peak_mko /peak_tko 
residues(22) = (simData_proTNF(1201,1)/max_val(1) - 0.2)/.04;
residues(23) = (simData_proTNF(1201,1)/simData_proTNF(1201,3) - 3.9)/.65;
residues(24) = (simData_proTNF(1201,1)/simData_proTNF(1201,2) - 14)/2.3; % wt_120 /
                                                          % tko_120
residues(25) = (simData_secTNF(301,3)/simData_secTNF(301,2) - 0.75)/.12;
% features of SecTNF
[max_val, max_ind] = max(simData_secTNF);
max_time = (max_ind-1)*0.1;

residues(26) = (simData_secTNF(601,3)/simData_secTNF(601,1) - 0.17)/.03;
residues(27) = (simData_secTNF(601,2)/simData_secTNF(601,3) - 1.4)/.54;
residues(28) = (simData_secTNF(1201,3)/simData_secTNF(1201,1) - .17)/.03;
residues(29) = (simData_secTNF(1201,2)/simData_secTNF(1201,3) - 1.9)/.76;
residues(30) = (simData_secTNF(1201,1)/simData_secTNF(601,1) - 1.2)/.16;


%
  inds = find(abs(residues)>1);
  residues(inds) = ones(size(inds))*999;


