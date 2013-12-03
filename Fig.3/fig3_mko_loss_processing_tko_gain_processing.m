% read data from R 
nascent_all = csvread('nascent2.csv',1,1);
nascent_exp = csvread('nascent.csv',1,0);
nascent_all = nascent_exp;
nascent_all(:,2:end) = nascent_exp(:,2:end)/k;
nascent_all(:,4) = nascent_all(:,4)/k;
nascent_all(:,6) = nascent_all(:,6)*2;

%kdeg_all = csvread('kdeg_all.csv',1,1);
yinit_all = csvread('yinit_all.csv',1,1);

times = 0:.1:120;%nascent_all(:,1);
[t,mRNA_wt]= ode15s(@ode,times,yinit_all(1),[],[],nascent_all(:,1:2), ...
                    kdeg_all(1));
[~,mRNA_mko]= ode15s(@ode,times,yinit_all(2),[],[],nascent_all(:,[1 ...
                    4]),kdeg_all(3));
[~,mRNA_tko]= ode15s(@ode,times,yinit_all(3),[],[],nascent_all(:,[1 ...
                    6]),kdeg_all(5));

csvwrite('mRNA_all_mko_loss_processing_tko_gain_processing.csv',[t mRNA_wt mRNA_mko mRNA_tko])
plot(t,mRNA_wt,'k')
hold on 
plot(t,mRNA_tko,'c')
plot(t,mRNA_mko,'Color', [0.5 0 0.5])
hold off
%function dydt = ode(t,x,options,nascent,kdeg)
%dydt= interp1(nascent(:,1),nascent(:,2),t) -kdeg*x ;
%end
