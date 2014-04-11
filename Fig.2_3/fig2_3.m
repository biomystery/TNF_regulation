% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
nascent_exp = csvread('../expdata/nascent.csv',1,0);
mRNA_exp = csvread('../expdata/mRNA.csv',1,0);

plot_flag = 0; 

pars = getParams(); % wt parameters
pars('pr_fold_mko') = 4.5;
pars('pr_fold_tko') = 1.1;
k_pr_all = [pars('k_pr') pars('k_pr')/(pars('pr_fold_mko')) pars('k_pr')/pars('pr_fold_tko')]; 
kdeg = [.02 .02 .07]; % wt, mko, tko 

yinit_all = zeros(2,3);
yinit_all(1,:) = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                  2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;
yinit_all(2,:) = yinit_all(1,:).*k_pr_all/kdeg(3);

times = 0:.1:120;%nascent_all(:,1);

% wt 
pars('kdeg_m') = kdeg(1);
[t,wt]= ode15s(@ode23,times,yinit_all(:,1),[],[],nfkb_exp(:,1:2), ...
               pars);
% mko 
pars('k_pr') = k_pr_all(2);
back_km = pars('Km_tr');
pars('Km_tr') = pars('Km_tr') *2; 
pars('kdeg_m') = kdeg(2);
[~,mko]= ode15s(@ode23,times,yinit_all(:,2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('k_pr') = k_pr_all(3);
pars('kdeg_m') = kdeg(3);
pars('Km_tr') = back_km; 
[~,tko]= ode15s(@ode23,times,yinit_all(:,3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);

figure
subplot 221
plot(t,wt(:,1),'k')
hold on 
plot(t,tko(:,1),'c')
plot(t,mko(:,1),'Color', [0.5 0 0.5])
legend('wt','tko','mko')
hold off
xlim([0 120])
title('nascent sim')

subplot 222
errorbar(nascent_exp(:,1),nascent_exp(:,2),nascent_exp(:,3),'k')
hold on 
errorbar(nascent_exp(:,1),nascent_exp(:,4),nascent_exp(:,5),'Color',[0.5 ...
                    0 0.5])
errorbar(nascent_exp(:,1),nascent_exp(:,6),nascent_exp(:,3),'c')
hold off
xlim([0 120])
title('nascent exp')

subplot 223 
plot(t,wt(:,2),'k')
hold on 
plot(t,tko(:,2),'c')
plot(t,mko(:,2),'Color', [0.5 0 0.5])
hold off
xlim([0 120])
title('mRNA sim')

subplot 224
errorbar(mRNA_exp(:,1),mRNA_exp(:,2),mRNA_exp(:,3),'k')
hold on 
errorbar(mRNA_exp(:,1),mRNA_exp(:,4),mRNA_exp(:,5),'Color',[0.5 ...
                    0 0.5])
errorbar(mRNA_exp(:,1),mRNA_exp(:,6),mRNA_exp(:,3),'c')
hold off
xlim([0 120])
title('mRNA exp')


%% save 
saveas(gca,'fig2_3.pdf')
close; 
!open fig2_3.pdf