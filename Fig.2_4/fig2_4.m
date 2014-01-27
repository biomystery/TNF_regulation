%% init 
% read expdata
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
nascent_exp = csvread('../expdata/nascent.csv',1,0);
mRNA_exp = csvread('../expdata/mRNA.csv',1,0);
sec_exp = csvread('../expdata/TNF_secrection.csv',1,0);
pro_exp = csvread('../expdata/combineTNF.csv',1,0);

% read pars
pars = getParams(); % wt parameters

% init condtions:
k_pr_all = [pars('k_pr') pars('k_pr')/(pars('pr_fold')*3) pars('k_pr')/pars('pr_fold')]; 
kdeg = [.02 .02 .07]; % wt, mko, tko 
k_tls = [pars('k_tl'),pars('k_tl'),pars('k_tl')/pars('tl_fold')]; 

yinit_all = zeros(2,3);
yinit_all(1,:) = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                  2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;
yinit_all(2,:) = yinit_all(1,:).*k_pr_all/kdeg(3);
yinit_all(3,:) = yinit_all(2,:).*k_tls/pars('kdeg_p');

% time
times = 0:.1:120;%nascent_all(:,1);

%% simulating 
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

proTNF = [wt(:,end) mko(:,end) tko(:,end)];

% secretion 
pars = getParams(); % set as default
k_secs = [pars('k_sec'),pars('k_sec'),pars('k_sec')/ ...
          pars('sec_fold')];

yinit_all(3,:) = yinit_all(2,:).*k_tls./(pars('kdeg_p')+k_secs); % init 

% wt 
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


secTNF = cumsum([wt(:,end)*k_secs(1) mko(:,end)*k_secs(2) tko(:,end)*k_secs(3)]);

%% plot
figure
subplot 241
plot(t,wt(:,1),'k')
hold on 
plot(t,tko(:,1),'c')
plot(t,mko(:,1),'Color', [0.5 0 0.5])

hold off
xlim([0 120])
title('nascent sim')

subplot 242
errorbar(nascent_exp(:,1),nascent_exp(:,2),nascent_exp(:,3),'k')
hold on 
errorbar(nascent_exp(:,1),nascent_exp(:,4),nascent_exp(:,5),'Color',[0.5 ...
                    0 0.5])
errorbar(nascent_exp(:,1),nascent_exp(:,6),nascent_exp(:,3),'c')
hold off
xlim([0 120])
title('nascent exp')

subplot 243 
plot(t,wt(:,2),'k')
hold on 
plot(t,tko(:,2),'c')
plot(t,mko(:,2),'Color', [0.5 0 0.5])
hold off
xlim([0 120])
title('mRNA sim')

subplot 244
errorbar(mRNA_exp(:,1),mRNA_exp(:,2),mRNA_exp(:,3),'k')
hold on 
errorbar(mRNA_exp(:,1),mRNA_exp(:,4),mRNA_exp(:,5),'Color',[0.5 ...
                    0 0.5])
errorbar(mRNA_exp(:,1),mRNA_exp(:,6),mRNA_exp(:,3),'c')
hold off
xlim([0 120])
title('mRNA exp')

% proTNF 
subplot 245
plot(t,proTNF(:,1),'k')
hold on 
plot(t,proTNF(:,3),'c') % tko 
plot(t,proTNF(:,2),'Color', [0.5 0 0.5]) %mko 
hold off
xlim([0 120])
title('proTNF sim')

subplot 246
plot(pro_exp(:,1),pro_exp(:,2),'.-k')
hold on 
plot(pro_exp(:,1),pro_exp(:,3),'.-','Color',[0.5 ...
                    0 0.5])
plot(pro_exp(:,1),pro_exp(:,4),'.-c')
hold off
xlim([0 120])
title('proTNF exp')

% sec TNF
subplot 247
plot(t,secTNF(:,1),'k')
hold on 
plot(t,secTNF(:,3),'c') %tko
plot(t,secTNF(:,2),'Color', [0.5 0 0.5]) % mko
hold off
xlim([0 120])
title('secTNF sim')
legend('wt','tko','mko')
subplot 248
errorbar(sec_exp(:,1),sec_exp(:,2),sec_exp(:,3),'k')
hold on 
errorbar(sec_exp(:,1),sec_exp(:,4),sec_exp(:,5),'Color',[0.5 ...
                    0 0.5])
errorbar(sec_exp(:,1),sec_exp(:,6),sec_exp(:,3),'c')
hold off
xlim([0 120])
title('secTNF exp')

%% save 
saveas(gca,'fig2_4.pdf')
close; 
!open fig2_4.pdf