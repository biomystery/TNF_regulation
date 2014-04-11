
% read data from R 
addpath('../src/')
mRNA_exp = csvread('../expdata/mRNA.csv',1,0);
expData = csvread('../expdata/proTNF.csv',1,0);
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

plot_flag = 1; 
N = 19;
degp_rate =[linspace(0.01,.09,9) linspace(1,10,10)];
sec_rate =[linspace(0.01,.09,9) linspace(1,10,10)];

for i =1:N % different pr vs same pr. 
    for j = 1:N
        input_pars(1) = degp_rate(i);
        input_pars(2) = sec_rate(j);
        residues = calScore_wt(input_pars,mRNA_exp,expData,0);
        score(i,j) = sum(residues.^2);
    end
    disp(i)
end
score = sqrt(score/numel(residues));
minimal_score = min(score(:));

figure('units','inch','position',[12 6 6 3])
subplot 121
if plot_flag
    plot(degp_rate,score,'k','linewidth',1.5)
    xlabel('kdeg_{p}');ylabel('RMSD');
    ylim([.2 .5])
    set(gca,'xscale','log','xgrid','on','ygrid','on')
end 

% find minimal 
minimal_score = min(score);
j = find(score == minimal_score);
subplot 122; calScore_wt(degp_rate(j),mRNA_exp,expData,1);
xlim([0 121]);ylim([0 1.4])
xlabel('Time (min)');ylabel('proTNF (a.u.)')

save ./simData/fig4s_proTNF_wt_partialTACE.mat 