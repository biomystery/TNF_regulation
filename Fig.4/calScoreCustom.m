function residues = calScoreCustom(input_pars,mRNA_all)

% params
plot_flag = 1;
pars = getParams(); % wt parameters

% input pars
pars('kdeg_p') = input_pars(1);
pars('k_sec') = input_pars(2);

%% simulations
% proTNF
k_tls = [pars('k_tl'),pars('k_tl'),pars('k_tl')/pars('tl_fold')]; 
yinit_all = mRNA_all(1,[2,4,6]).*k_tls/pars('kdeg_p'); % init 

times = 0:.1:120;%max(mRNA_all(:,1));%mRNA_all(:,1);

[t,proTNF_wt]= ode15s(@ode4,times,yinit_all(1),[],[],mRNA_all(:,1:2), ...
                    pars);
[~,proTNF_mko]= ode15s(@ode4,times,yinit_all(2),[],[],mRNA_all(:,[1 ...
                    4]),pars);

pars('k_tl') = k_tls(3); % set as tko
[~,proTNF_tko]= ode15s(@ode4,times,yinit_all(3),[],[],mRNA_all(:,[1 ...
                    6]),pars);
pars('k_tl') = k_tls(1); % set as default
simData_proTNF = [proTNF_wt proTNF_mko proTNF_tko];

% secTNF

k_secs = [pars('k_sec'),pars('k_sec'),pars('k_sec')/pars('sec_fold')];
yinit_all = mRNA_all(1,[2,4,6]).*k_tls./(pars('kdeg_p')+k_secs); % init 


[t,TNF_wt]= ode15s(@ode4_1,times,yinit_all(1),[],[],mRNA_all(:,1:2), ...
                    pars);
[~,TNF_mko]= ode15s(@ode4_1,times,yinit_all(2),[],[],mRNA_all(:,[1 ...
                    4]),pars);

pars('k_tl') = k_tls(3); % set as tko
pars('k_sec') = k_secs(3); % set as tko
[~,TNF_tko]= ode15s(@ode4_1,times,yinit_all(3),[],[],mRNA_all(:,[1 ...
                    6]),pars);
pars('k_tl') = k_tls(1); % set as default
pars('k_sec') = k_secs(1); % set as default

simData_secTNF = [cumsum([TNF_wt*k_secs(1) TNF_mko*k_secs(2) TNF_tko*k_secs(3)])];


% plot 
if plot_flag
    figure
    subplot(2,1,1)
    plot(times,simData_proTNF)
    xlim([0 120])
    subplot(2,1,2)
    plot(times,simData_secTNF)
    xlim([0 120])
end

%% calculate score 
% features of proTNF
residues = zeros(1,14); % 

[max_val, max_ind] = max(simData_proTNF);
max_time = (max_ind-1)*0.1;

residues(1) = (max_time(1) - 60)/10 ; % wt peak time 
residues(2) = (max_time(2) - 60)/10 ; % mko peak time 
residues(3) = (max_time(3) - 60)/10 ; % tko peak time 
residues(4) = (max_val(3) / max_val(1) - 0.34)/0.06; % peak_tko / peak_wt 
residues(5) = (max_val(2) / max_val(3) - 0.82)/0.14; % peak_mko /peak_tko 
residues(6) = (simData_proTNF(1201,1)/max_val(1) - 0.2)/.04;
residues(7) = (simData_proTNF(1201,1)/simData_proTNF(1201,3) - 3.9)/.65;
residues(8) = (simData_proTNF(1201,1)/simData_proTNF(1201,2) - 14)/2.3; % wt_120 /
                                                          % tko_120
residues(9) = (simData_secTNF(301,3)/simData_secTNF(301,2) - 0.75)/.12;
% features of SecTNF
[max_val, max_ind] = max(simData_secTNF);
max_time = (max_ind-1)*0.1;

residues(10) = (simData_secTNF(601,3)/simData_secTNF(601,1) - 0.17)/.03;
residues(11) = (simData_secTNF(601,2)/simData_secTNF(601,3) - 1.4)/.54;
residues(12) = (simData_secTNF(1201,3)/simData_secTNF(1201,1) - .17)/.03;
residues(13) = (simData_secTNF(1201,2)/simData_secTNF(1201,3) - 1.9)/.76;
residues(14) = (simData_secTNF(1201,1)/simData_secTNF(601,1) - 1.2)/.16;



