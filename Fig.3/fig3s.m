% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
expData = csvread('../expdata/mRNA.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 

% plot options 
plot_flag = 1; 

N = 10;
pr_fold = linspace(1,10,N);
rmsd = zeros(N,N);
nrmsd = zeros(N,N);

for i = 1:N
    for j = 1:N

        pars = getParams(); % wt parameters
        k_pr_all = [pars('k_pr') pars('k_pr')/pr_fold(i) pars('k_pr')/pr_fold(j)]; 
        kdeg_m = [.02 .02 .07]; % wt, mko, tko 
        k_pr = pars('k_pr');
        kdeg_m_mko = [0.07, 0.02];

        % initial conditions
        yinit = nascent_exp(1,2:end) .* k_pr_all/kdeg_m(3);
        times = 0:.1:120;%nascent_all(:,1);

        % wt 
        [t,wt]= ode15s(@ode3,times,yinit(:,1),[],[],nascent_exp(:,1:2), ...
                       pars);
        pars('k_pr') = k_pr_all(1);
        pars('kdeg_m') = kdeg_m(1);

        % mko 
        pars('k_pr') = k_pr_all(2);
        pars('kdeg_m') = kdeg_m(2);
        [~,mko]= ode15s(@ode3,times,yinit(:,2),[],[],nascent_exp(:,[1 ...
                            3]),pars);

        % tko 
        pars('k_pr') = k_pr_all(3);
        pars('kdeg_m') = kdeg_m(3);
        [~,tko]= ode15s(@ode3,times,yinit(:,3),[],[],nascent_exp(:,[1 ...
                            4]),pars);
        
        % sim data 
        simData = [wt mko tko];
        simData = simData((expData(1:4,1))*10+1,:);

        % calculate score 
        [rmsd(i,j), nrmsd(i,j)] = calNRMSD(simData,expData(1:4,2:end));
        disp(((i-1)*N+j)/(N*N))
        
    end
end

%% R = RMSD / RMSD_lin
expData_mean = expData(:,[2,4,6]);
expData_mean = expData_mean(:);
rmsd_lin = sqrt(sum((expData_mean-mean(expData_mean)).^2)/ ...
                numel(expData_mean));



% write the data into the file 
%csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])
if plot_flag
figure
imagesc(pr_fold,pr_fold,rmsd/rmsd_lin);colorbar
xlabel('pr_{fold} (tko)')
ylabel('pr_{fold} (mko)')
end