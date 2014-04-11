% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData =  csvread('../expdata/nascent.csv',1,0);
plot_flag = 1; 
pars = getParams(); % wt parameters
pr_fold = [pars('pr_fold') 1];
k_pr = pars('k_pr');


N = 30;
pr_fold = linspace(0.1,3,N);
rmsd = zeros(N,N);
nrmsd = zeros(N,N);

for i =1:N
    for j = 1:N

        k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold(i) pars('k_pr')/pr_fold(j)]; 
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
    
        % get date
        simData = [nascent_wt nascent_mko nascent_tko];
        simData = simData((expData(:,1))*10+1,:);

        % calculate score 
        [rmsd(i,j), nrmsd(i,j)] = calNRMSD(simData,expData(:,2:end));
        disp(((i-1)*N+j)/(N*N))
    end
end

save ./simData/fig2s.mat 

% plot 
if plot_flag
    imagesc(pr_fold,pr_fold,nrmsd)
end 

% save data 
csvwrite('./simData/nrmsd.csv',nrmsd)

% run R
!R CMD BATCH fig2s.R
!rm *.Ro* 
!rm *.Rh*
!rm *.RD*