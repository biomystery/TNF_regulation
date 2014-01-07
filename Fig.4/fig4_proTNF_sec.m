% read data from R 
mRNA_all = csvread('../expdata/mRNA.csv',1,0);
mRNA_all(:,2:end) = mRNA_all(:,2:end)/9/15*8.4; % tl_rate
kdeg = 5.8e-3*10; % from Werner et al. 2008
tl_regulator  = 1.5; %  fold less in tko or 1.5
mRNA_all(:,[2,4,6]) = mRNA_all(:,[2,4,6])/tl_regulator; % only TKO has slower
                                                % translation rate 

yinit_all = mRNA_all(1,[2,4,6])/kdeg;

times = 0:.1:max(mRNA_all(:,1));%mRNA_all(:,1);
[t,proTNF_wt]= ode15s(@ode,times,yinit_all(1),[],[],mRNA_all(:,1:2), ...
                    kdeg);
[~,proTNF_mko]= ode15s(@ode,times,yinit_all(2),[],[],mRNA_all(:,[1 ...
                    4]),kdeg);
[~,proTNF_tko]= ode15s(@ode,times,yinit_all(3),[],[],mRNA_all(:,[1 ...
                    6]),kdeg);

csvwrite('proTNF_all_sec.csv',[t proTNF_wt proTNF_mko proTNF_tko])


% secretion 

ksec = kdeg;
sec_regulator = 5; % 2.5
                     %ksec = kdeg/sec_regulator;

yinit_all = mRNA_all(1,[2,4,6])/(kdeg + ksec);


[t,TNF_wt]= ode15s(@ode2,times,yinit_all(1),[],[],mRNA_all(:,1:2), ...
                    kdeg,ksec);
[~,TNF_mko]= ode15s(@ode2,times,yinit_all(2),[],[],mRNA_all(:,[1 ...
                    4]),kdeg,ksec);
[~,TNF_tko]= ode15s(@ode2,times,yinit_all(3),[],[],mRNA_all(:,[1 ...
                    6]),kdeg,ksec/sec_regulator);

csvwrite('sec_all_sec.csv',[t cumsum([TNF_wt*ksec TNF_mko*ksec TNF_tko*ksec/sec_regulator])]);


%end
%!R CMD BATCH Fig4_single.R
