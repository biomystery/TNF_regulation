function residues = calScoreCustom23(input_pars)

% expdata
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);

% params
plot_flag = 1;


pars = getParams(); % wt parameters
input_pars = 10.^input_pars;
pars('V_tr') = 1;%input_pars(6); 
pars('Km_tr') = input_pars(1);
pars('n') = input_pars(2); 

kdeg = [.02 .02 .07]; % wt, mko, tko 

% init     
yinit_all = zeros(1,3);
yinit_all(1,:) = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                  2:end).^pars('n')+pars('Km_tr')^pars('n'))./kdeg(3);
% time 
times = 0:.1:120;%nascent_all(:,1);

%% simulations
% wt
[t,wt]= ode15s(@ode23_new,times,yinit_all(:,1),[],[],nfkb_exp(:,1:2), ...
                       pars);

% mko 
pars('kdeg_m') = kdeg(2);
[~,mko]= ode15s(@ode23_new,times,yinit_all(:,2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('kdeg_m') = kdeg(3);
[~,tko]= ode15s(@ode23_new,times,yinit_all(:,3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);

% get date
simData_mRNA = [wt(:,1) mko(:,1) tko(:,1)];


if plot_flag
    subplot(2,2,1)
    plot(times,simData_mRNA)
    hold off
end


%% calculate score 
% residue1, peak time 
residues = zeros(1,16); % 8 features. 


% features of mRNA profile
[max_val, max_ind] = max(simData_mRNA);
max_time = (max_ind-1)*0.1;

residues(1) = (max_time(1) - 60)/10 ; % wt peak time 
residues(2) = (max_time(2) - 60)/20 ; % mko peak time 
residues(3) = (max_time(3) - 60)/10 ; % tko peak time 

residues(4) = (max_val(3) / max_val(1) - 0.73)/0.05; % peak_tko / peak_wt 
residues(5) = (max_val(2) / max_val(3) - 0.5)/0.06; % peak_mko /peak_tko 
residues(6) = (simData_mRNA(1201,1)/simData_mRNA(1201,2) - 1.98)/.23;
residues(7) = (simData_mRNA(1201,1)/simData_mRNA(1201,3) - 1.17)/.14;
