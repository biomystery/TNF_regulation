% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);

plot_flag = 0; 

pars = getParams(); % wt parameters

k_pr_all = [pars('V_pr') pars('V_pr')/1.5 pars('V_pr')/1.5]; 
kdeg = [.02 .02 .07]; % wt, mko, tko 

yinit_all = zeros(2,3);
yinit_all(1,:) = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                  2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;
yinit_all(2,:) = yinit_all(1,:).*k_pr_all/kdeg(3);

times = 0:.1:120;%nascent_all(:,1);

% wt 
[t,wt]= ode15s(@ode23_new,times,yinit_all(:,1),[],[],nfkb_exp(:,1:2), ...
               pars);

% mko 
pars('V_pr') = k_pr_all(2);
pars('kdeg_m') = kdeg(2);
[~,mko]= ode15s(@ode23_new,times,yinit_all(:,2),[],[],nfkb_exp(:,[1 ...
                    3]),pars);

% tko 
pars('V_pr') = k_pr_all(3);
pars('kdeg_m') = kdeg(3);

[~,tko]= ode15s(@ode23_new,times,yinit_all(:,3),[],[],nfkb_exp(:,[1 ...
                    4]),pars);

figure
subplot 211
plot(t,wt(:,1),'k')
hold on 
plot(t,tko(:,1),'c')
plot(t,mko(:,1),'Color', [0.5 0 0.5])
hold off
subplot 212
plot(t,wt(:,2),'k')
hold on 
plot(t,tko(:,2),'c')
plot(t,mko(:,2),'Color', [0.5 0 0.5])
hold off

saveas(gca,'fig2_3.pdf')

%function dydt = ode(t,x,options,nascent,kdeg)
%dydt= interp1(nascent(:,1),nascent(:,2),t) -kdeg*x ;
%end
