set(0,'DefaultAxesColorOrder',[54 100 139; 135 206 255;79 148 205]/255)
% read data from R 
addpath('../src/')
mRNA_exp = csvread('../expdata/mRNA.csv',1,0);
mRNA_exp(:,[3,5,7]) = [];
expData = csvread('../expdata/combineTNF.csv',1,0);
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

expData2 = csvread('../expdata/TNF_secrection.csv',1,0);

expData2(:,2:end) = expData2(:,2:end)/max(expData2(:,2));  % normalized expData
expData2(5:6,:) = [];

plot_flag = 1; 
N = 100;

degp_rate =linspace(.01,1,100);



for i =1:N % different pr vs same pr. 

    input_pars(1) = degp_rate(i);

    residues = calScore(input_pars,mRNA_exp,expData,0);
    score(i) = sum(residues.^2);
    disp(i)
end
score = sqrt(score/numel(residues));


figure('units','inch','position',[12 6 6 3])
subplot 121
if plot_flag
    plot(degp_rate,score,'linewidth',1.5)
    xlabel('kdeg_{p}(min^{-1})');ylabel('RMSD');
end 

% find minimal 
minimal_score = min(score(:));
i= find((score) == minimal_score);

subplot 122; calScore(degp_rate(i(1)),mRNA_exp,expData,1);
xlim([0 121]);ylim([0 1.4])
xlabel('Time (min)');ylabel('proTNF (a.u.)')


text(20,1.3,strcat('RMSD=',num2str(minimal_score)))
text(20,1.15,strcat('kdeg_p=',num2str(degp_rate(i(1)))))

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

print -dpdf -painters fig4s_all_proTNF.pdf
saveas(gca,'fig4s_all_proTNF.fig')
close all; !open fig4s_all_proTNF.pdf

save ./simData/fig4s_all_proTNF.mat 