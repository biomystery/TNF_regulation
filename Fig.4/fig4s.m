% read data from R 
mRNA_all = csvread('../expdata/mRNA.csv',1,0);
sec_exp = csvread('../expdata/TNF_secrection.csv',1,0)
pro_exp = csvread('../expdata/combineTNF.csv',1,0)
        pro_exp = [pro_exp(:,1:2) zeros(5,1) pro_exp(:,3) zeros(5,1) pro_exp(:,4) ...
                   zeros(5,1)];

%% pro TNF
N = 20;
M = 20;
tl_fold = linspace(0.5,10,N);
ktls = linspace(0.05,1,M);
mRNA_back = mRNA_all; 
rmsd_pro = zeros(N,M);
nrmsd_pro = zeros(N,M);
%
for i = 1:N
    for j = 1:M
        
        %ktl = 1/9/15*8.4;
        %tl_regulator  = 1.5; %  fold less in tko or 1.5
        ktl = ktls(i);
        tl_regulator = tl_fold(j);
        mRNA_all = mRNA_back;

        mRNA_all(:,2:end) = mRNA_all(:,2:end)*ktl;
        kdeg = 5.8e-3*10; % from Werner et al. 2008


        mRNA_all(:,[6]) = mRNA_all(:,[6])/tl_regulator;

        yinit_all = mRNA_all(1,[2,4,6])/kdeg;

        times = 0:.1:max(mRNA_all(:,1));%mRNA_all(:,1);
        [t,proTNF_wt]= ode15s(@ode,times,yinit_all(1),[],[],mRNA_all(:,1:2), ...
                              kdeg);
        [~,proTNF_mko]= ode15s(@ode,times,yinit_all(2),[],[],mRNA_all(:,[1 ...
                            4]),kdeg);
        [~,proTNF_tko]= ode15s(@ode,times,yinit_all(3),[],[],mRNA_all(:,[1 ...
                            6]),kdeg);

        pro_sim = [proTNF_wt proTNF_mko proTNF_tko];
        %csvwrite('proTNF_all.csv',[t pro_sim])


        pro_sim = pro_sim(pro_exp(:,1)*10+1,:);
        [rmsd_pro(i,j), nrmsd_pro(i,j)] = calNRMSD(pro_sim,pro_exp(:,2:end));
        
    end
end

%% secretion 

N = 20;
M = 20;
sec_fold = linspace(0.5,10,N);
ksecs = linspace(0.05,1,M);
mRNA_back = mRNA_all; 
rmsd_sec = zeros(N,M);
nrmsd_sec = zeros(N,M);
%
mRNA_all = mRNA_back;
for i = 1:N
    for j = 1:M

        %sec_regulator = 5; % 2.5
        %ksec =  kdeg;
        ksec = ksecs(i);
        sec_regulator = sec_fold(j);
        yinit_all = mRNA_all(1,[2,4,6])/(kdeg + ksec);


        [t,TNF_wt]= ode15s(@ode2,times,yinit_all(1),[],[],mRNA_all(:,1:2), ...
                           kdeg,ksec);
        [~,TNF_mko]= ode15s(@ode2,times,yinit_all(2),[],[],mRNA_all(:,[1 ...
                            4]),kdeg,ksec);
        [~,TNF_tko]= ode15s(@ode2,times,yinit_all(3),[],[],mRNA_all(:,[1 ...
                            6]),kdeg,ksec/sec_regulator);

        sec_sim =cumsum([TNF_wt*ksec TNF_mko*ksec TNF_tko*ksec/ ...
                         sec_regulator]);


        sec_sim = sec_sim(sec_exp(1:4,1)*10+1,:);
        [rmsd_sec(i,j), nrmsd_sec(i,j)] = calNRMSD(sec_sim,sec_exp(1:4,2: ...
                                                          end));
        
    end
end 

%%
figure
imagesc(ktls,tl_fold,log10(rmsd_pro));colorbar
xlabel('k_{tl}')
ylabel('tl fold reduction (tko)')

figure 
imagesc(ksecs,sec_fold,log10(rmsd_sec));colorbar
xlabel('k_{sec}')
ylabel('sec fold reduction (tko)')
%end
%!R CMD BATCH Fig4_single.R
