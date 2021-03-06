%% figs1. The figure for same and different processing rate 

% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
plot_flag = 0; 
pars = getParams(); % wt parameters
pr_fold = [1.5 4.5]; 
k_pr = pars('k_pr');


% set up the file name to same the simulate data 
filenames = {'./simData/fig2s2.csv','./simData/fig2s2b.csv'};


for i =1:2 % different pr vs same pr. 

    k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold(2) ...
                pars('k_pr')/pr_fold(1)];
    
    yinit_all = pars('V_tr')* nfkb_exp(1,2:end).^pars('n')./(nfkb_exp(1, ...
                                                      2:end).^pars('n')+pars('Km_tr')^pars('n'))./k_pr_all;

    times = 0:.1:120;%nascent_all(:,1);

    % wt
    [t,nascent_wt]= ode15s(@ode2,times,yinit_all(1),[],[],nfkb_exp(:,1:2), ...
                           pars);

    % mko 
    pars('k_pr') = k_pr_all(2);

    if i == 2
        pars('Km_tr') = pars('Km_tr')*2;
    end
    
    [~,nascent_mko]= ode15s(@ode2,times,yinit_all(2),[],[],nfkb_exp(:,[1 ...
                        3]),pars);
    if i == 2
        pars('Km_tr') = pars('Km_tr')/2;
    end

    % tko 
    pars('k_pr') = k_pr_all(3);
    [~,nascent_tko]= ode15s(@ode2,times,yinit_all(3),[],[],nfkb_exp(:,[1 ...
                        4]),pars);

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
!R CMD BATCH fig2s2.R
!rm *.Ro* 
