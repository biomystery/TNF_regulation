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
pars('kdeg_p') = log(2)/10;
pars('k_sec') =  log(2)/10;
% init condtions:
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
%figure()
figure('Position', [100, 100, 700, 700])
subplot 221
[ax,h1,h2] = plotyy(t,wt(:,1),nascent_exp(:,1),nascent_exp(:,2),'plot');
set(h1,'Color','k')
set(h2,'Visible','off')
set(gcf, 'CurrentAxes', ax(1));
hold on 
plot(t,tko(:,1),'c')
plot(t,mko(:,1),'Color', [0.5 0 0.5])
xlim([0 120])
ylim([0 3])

set(gcf, 'CurrentAxes', ax(2));
ylim([0 33])
hold on 
errorbar(nascent_exp(:,1),nascent_exp(:,2),nascent_exp(:,3),'ko');
errorbar(nascent_exp(:,1),nascent_exp(:,4),nascent_exp(:,5),'o','Color',[0.5 ...
                    0 0.5])
errorbar(nascent_exp(:,1),nascent_exp(:,6),nascent_exp(:,7),'co')
xlim([0 120])

title('nascent')
set(ax(1),'xtick',0:30:120,'ytick',0:3)
set(ax(2),'xtick',0:30:120,'ytick',0:10:30)


subplot 222 
[ax,h1,h2] = plotyy(t,wt(:,2),mRNA_exp(:,1),mRNA_exp(:,2),'plot');
set(h1,'Color','k')
set(h2,'Visible','off')
set(gcf, 'CurrentAxes', ax(1));
hold on 
plot(t,tko(:,2),'c')
plot(t,mko(:,2),'Color', [0.5 0 0.5])
xlim([0 120])
ylim([0 20])

set(gcf, 'CurrentAxes', ax(2));
ylim([0 200])
hold on 
errorbar(mRNA_exp(:,1),mRNA_exp(:,2),mRNA_exp(:,3),'ko');
errorbar(mRNA_exp(:,1),mRNA_exp(:,4),mRNA_exp(:,5),'o','Color',[0.5 ...
                    0 0.5])
errorbar(mRNA_exp(:,1),mRNA_exp(:,6),mRNA_exp(:,7),'co')
xlim([0 120])

title('mRNA')
set(ax(1),'xtick',0:30:120,'ytick',0:5:20)
set(ax(2),'xtick',0:30:120,'ytick',0:40:200)


% proTNF 
subplot 223

[ax,h1,h2] = plotyy(t,proTNF(:,1),pro_exp(:,1),pro_exp(:,2),'plot');
set(h1,'Color','k')
set(h2,'LineStyle','none','Marker','o','Color','k')
set(gcf, 'CurrentAxes', ax(1));
hold on 
plot(t,proTNF(:,3),'c')
plot(t,proTNF(:,2),'Color', [0.5 0 0.5])
xlim([0 120])


set(gcf, 'CurrentAxes', ax(2));
ylim([0 120])
hold on 
plot(pro_exp(:,1),pro_exp(:,3),'o','Color',[0.5 ...
                    0 0.5])
plot(pro_exp(:,1),pro_exp(:,4),'oc')

xlim([0 120])
title('proTNF')

set(ax(1),'xtick',0:30:120)%,'ytick',0:.5:2)
set(ax(2),'xtick',0:30:120,'ytick',0:30:120)


subplot 224

[ax,h1,h2] = plotyy(t,secTNF(:,1),sec_exp(:,1),sec_exp(:,2),'plot');
set(h1,'Color','k')
set(h2,'Visible','off')

set(gcf, 'CurrentAxes', ax(1));
hold on 
plot(t,secTNF(:,3),'c')
plot(t,secTNF(:,2),'Color', [0.5 0 0.5])
xlim([0 120])

legend('wt','tko','mko','Location','best')
legend boxoff 

set(gcf, 'CurrentAxes', ax(2));

hold on 
errorbar(sec_exp(:,1),sec_exp(:,2),sec_exp(:,3),'ko')
errorbar(sec_exp(:,1),sec_exp(:,4),sec_exp(:,5),'o','Color',[0.5 ...
                    0 0.5])
errorbar(sec_exp(:,1),sec_exp(:,6),sec_exp(:,3),'co')

xlim([0 120])
title('secTNF')

set(ax(1),'xtick',0:30:120)%,'ytick',0:.5:2)
set(ax(2),'xtick',0:30:120)%,'ytick',0:30:120)

%% 
saveas(gca,'fig2_4b.pdf')
saveas(gca,'fig2_4b.fig')
close; 
!open fig2_4b.pdf