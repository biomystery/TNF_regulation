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
N = 10;
M = 10;

degp_rate =linspace(.1,1,10);
sec_rate =linspace(1,10,10);

score = zeros(N,M);


for i =1:N % different pr vs same pr. 
    for j = 1:M
    input_pars(1) = degp_rate(i);
    input_pars(2) = sec_rate(j);
    residues = calScoreb2(input_pars,mRNA_exp,expData2,0);
    score(i,j) = sum(residues.^2);
    end
    disp(i)
end
score = sqrt(score/numel(residues));


figure('units','inch','position',[12 6 8 3])
subplot 121
if plot_flag
    imagesc(sec_rate,degp_rate,score)
    ylabel('kdeg_{p}(min^{-1})');xlabel('k_{sec}fold_{tko}');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 
minimal_score = min(score(:));
j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 


subplot 122;calScoreb2([degp_rate(i),sec_rate(j)],mRNA_exp,expData2,1);
xlim([0 121]);ylim([0 1.4])
xlabel('Time (min)');ylabel('secTNF (a.u.)')

text(20,1.3,strcat('RMSD=',num2str(minimal_score)))
text(20,1.15,strcat('kdeg_p=',num2str(degp_rate(i(1)))))
text(20,1,strcat('k_{sec}fold_{tko}=',num2str(sec_rate(j(1)))))
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

print -dpdf -painters fig4s_all_secTNF2.pdf
saveas(gca,'fig4s_all_secTNF2.fig')
close all; !open fig4s_all_secTNF2.pdf

save ./simData/fig4s_all_secTNF2.mat 