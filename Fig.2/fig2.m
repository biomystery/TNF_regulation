% read data from R 
nfkb_exp = csvread('nfkb_input.csv',1,0);
nfkb_exp(:,2)=nfkb_exp(:,2); % wt less induction
k = [.15 .1 .1]*4; 
yinit_all = nfkb_exp(1,2:end).^3./(nfkb_exp(1,2:end).^3+.5^3)./k;

times = 0:.1:120;%nascent_all(:,1);

[t,nascent_wt]= ode15s(@ode,times,yinit_all(1),[],[],nfkb_exp(:,1:2), ...
                    k(1));
[~,nascent_mko]= ode15s(@ode,times,yinit_all(2),[],[],nfkb_exp(:,[1 ...
                    3]),k(2));
[~,nascent_tko]= ode15s(@ode,times,yinit_all(3),[],[],nfkb_exp(:,[1 ...
                    4]),k(3));

csvwrite('nascent_mko_1.5_fold_less_process.csv',[t nascent_wt nascent_mko ...
                    nascent_tko])
figure
plot(t,nascent_wt,'k')
hold on 
plot(t,nascent_tko,'c')
plot(t,nascent_mko,'Color', [0.5 0 0.5])
hold off
%function dydt = ode(t,x,options,nascent,kdeg)
%dydt= interp1(nascent(:,1),nascent(:,2),t) -kdeg*x ;
%end
