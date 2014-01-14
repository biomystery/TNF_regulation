% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);

plot_flag = 0; 

pars = getParams(); % wt parameters
pr_fold = [pars('pr_fold') 1];
k_pr = pars('k_pr');
filenames = {'./simData/different_pr.csv','./simData/same_pr.csv'};
for i =1:2
    pars('pr_fold') = pr_fold(i);
    k_pr_all = [pars('k_pr') pars('k_pr')/pars('pr_fold') pars('k_pr')/pars('pr_fold')]; 
    yinit_all = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                      2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;

    times = 0:.1:120;%nascent_all(:,1);

    % wt
    [t,nascent_wt]= ode15s(@ode2,times,yinit_all(1),[],[],nfkb_exp(:,1:2), ...
                           pars);

    % mko 
    pars('k_pr') = k_pr_all(2);
    [~,nascent_mko]= ode15s(@ode2,times,yinit_all(2),[],[],nfkb_exp(:,[1 ...
                        3]),pars);

    % tko 
    pars('k_pr') = k_pr_all(3);
    [~,nascent_tko]= ode15s(@ode2,times,yinit_all(3),[],[],nfkb_exp(:,[1 ...
                        4]),pars);

    pars('k_pr') = k_pr;
    % save the date
    csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])

    % draw the plot 
    if plot_flag
        figure
        plot(t,nascent_wt,'k')
        hold on 
        plot(t,nascent_tko,'c')
        plot(t,nascent_mko,'Color', [0.5 0 0.5])
        hold off
    end
end

% run R
!R CMD BATCH fig2.R
!rm *.Ro* 
!rm *.Rh*
!rm *.RD*