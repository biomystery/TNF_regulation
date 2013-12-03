% read data from R 
nfkb_exp = csvread('nfkb_input.csv',1,0);
k = [.15 .1 .1]*4; % wt, mko, tko
kdeg = [.02 .02 .07]; % wt, mko, tko 
yinit_all = zeros(2,3);
yinit_all(1,:) = 27/6* nfkb_exp(1,2:end).^3./(nfkb_exp(1,2:end).^3+.5^3)./ ...
    k;
yinit_all(2,:) = yinit_all(1,:).*k./kdeg;

times = 0:.1:120;%nascent_all(:,1);

[t,wt]= ode15s(@ode23,times,yinit_all(:,1),[],[],nfkb_exp(:,1:2), ...
                    k(1),kdeg(1));
[~,mko]= ode15s(@ode23,times,yinit_all(:,2),[],[],nfkb_exp(:,[1 ...
                    3]),k(2),kdeg(2));
[~,tko]= ode15s(@ode23,times,yinit_all(:,3),[],[],nfkb_exp(:,[1 ...
                    4]),k(3),kdeg(3));

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


%function dydt = ode(t,x,options,nascent,kdeg)
%dydt= interp1(nascent(:,1),nascent(:,2),t) -kdeg*x ;
%end
