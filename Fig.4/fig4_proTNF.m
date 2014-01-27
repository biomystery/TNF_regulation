% read data from R 
mRNA_all = csvread('../expdata/mRNA.csv',1,0);
pars = getParams(); 


k_tls = [pars('k_tl'),pars('k_tl'),pars('k_tl')]/pars('tl_fold'); 
yinit_all = mRNA_all(1,[2,4,6]).*k_tls/pars('kdeg_p'); % init 


times = 0:.1:max(mRNA_all(:,1));%mRNA_all(:,1);
back_tl = pars('k_tl');
pars('k_tl') = k_tls(1);
[t,proTNF_wt]= ode15s(@ode4,times,yinit_all(1),[],[],mRNA_all(:,1:2), ...
                    pars);
[~,proTNF_mko]= ode15s(@ode4,times,yinit_all(2),[],[],mRNA_all(:,[1 ...
                    4]),pars);
[~,proTNF_tko]= ode15s(@ode4,times,yinit_all(3),[],[],mRNA_all(:,[1 ...
                    6]),pars);
%pars('k_tl') = back_tl;

csvwrite('./simData/proTNF_all.csv',[t proTNF_wt proTNF_mko proTNF_tko])


% secretion 
k_secs = [pars('k_sec'),pars('k_sec'),pars('k_sec')]/pars('sec_fold');
yinit_all = mRNA_all(1,[2,4,6]).*k_tls./(pars('kdeg_p')+k_secs); % init 

back_sec = pars('k_sec'); 
pars('k_sec') = k_secs(1);

[t,TNF_wt]= ode15s(@ode4_1,times,yinit_all(1),[],[],mRNA_all(:,1:2), ...
                    pars);
[~,TNF_mko]= ode15s(@ode4_1,times,yinit_all(2),[],[],mRNA_all(:,[1 ...
                    4]),pars);
[~,TNF_tko]= ode15s(@ode4_1,times,yinit_all(3),[],[],mRNA_all(:,[1 ...
                    6]),pars);

csvwrite('./simData/sec_all.csv',[t cumsum([TNF_wt*k_secs(1) TNF_mko*k_secs(2) TNF_tko*k_secs(3)])]);


%end
%!R CMD BATCH Fig4_single.R
