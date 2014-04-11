
% read data from R 
addpath('../src/')
mRNA_exp = csvread('../expdata/mRNA.csv',1,0);
mRNA_exp(:,[3,5,7]) = [];
expData = csvread('../expdata/combineTNF.csv',1,0);
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

plot_flag = 1; 
N = 19;
M = 19;

degp_rate =[linspace(0.01,.09,9) linspace(.1,1,10)];
sec_rate =[linspace(0.01,.09,9) linspace(.1,1,10)];

score = zeros(N,M);

for i =1:N % different pr vs same pr. 
    for j = 1:M
    input_pars(1) = degp_rate(i);
    input_pars(2) = sec_rate(j);
    residues = calScore(input_pars,mRNA_exp,expData,0);
    score(i,j) = sum(residues.^2);
    end
    disp(i)
end
score = sqrt(score/numel(residues));


figure('units','inch','position',[12 6 7 3])
subplot 121
if plot_flag
    imagesc((sec_rate),(degp_rate),score,[minimal_score,minimal_score+.1])
    ylabel('log10(kdeg_{p}(min^{-1}))');xlabel('k_{sec}(min^{-1})');
    h = colorbar;
    ylabel(h,'RMSD')
end 


% find minimal 
minimal_score = min(score(:));
j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 

subplot 122; calScore([degp_rate(i(1)),sec_rate(j(1))],mRNA_exp,expData,1);
xlim([0 121]);ylim([0 1.4])
xlabel('Time (min)');ylabel('proTNF (a.u.)')

text(20,1.3,strcat('RMSD=',num2str(minimal_score)))
text(20,1.15,strcat('kdeg_p=',num2str(degp_rate(i(1)))))
text(20,1,strcat('k_{sec}=',num2str(sec_rate(j(1)))))

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fig4s_proTNF_partialTace.pdf
saveas(gca,'fig4s_proTNF_partialTace.fig')
close all; !open fig4s_proTNF_partialTace.pdf

save ./simData/fig4s_proTNF_partialTace.mat 