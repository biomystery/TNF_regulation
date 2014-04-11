%% figs1. The figure for same and different processing rate 
set(0,'DefaultAxesColorOrder',[54 100 139; 135 206 255;79 148 205]/255)
% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
nascent_exp = csvread('../expdata/nascent.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 

expData1 = csvread('../expdata/nascent.csv',1,0);
expData1(:,2:end) = expData1(:,2:end)/max(expData1(:,2));  %
                                                           % normalized expData

expData2 = csvread('../expdata/mRNA.csv',1,0);
expData2 = expData2(1:4,:); 
expData2(1,3) = 0.01;
expData2(:,2:end) = expData2(:,2:end)/max(expData2(:,2));  % normalized expData

plot_flag = 1; 
pars = getParams(); % wt parameters

N = 41;
pr_rate = linspace(1,5,N);


for i =1:N % different pr vs same pr. 
    input_pars(1) = pr_rate(i);
    residues = calScores_wtmko(input_pars,nfkb_exp,expData1,0);
    score(i) = sum(residues.^2);
    disp(i)
end

score = sqrt(score/numel(residues));


minimal_score = min(score(:));
figure('units','inch','position',[12 6 6 3])
if plot_flag
    subplot 121
    plot(pr_rate,score,'k','linewidth',1.5)
    ylabel('RMSD');xlabel('k_{tr}fold_{mko}');

end 
xlim([1,5])
set(gca,'xgrid','on')

% find minimal 

i = find(score == minimal_score);

subplot 122
calScores_wtmko(pr_rate(i(1)),nfkb_exp,expData1,1);
text(50,1.8,strcat('RMSD=',num2str(minimal_score))
text(50,1.6,strcat('k_{tr}fold_{mko}=',num2str(pr_rate(i))))

xlim([0 121]);ylim([0 2])
xlabel('Time (min)');ylabel('nascent (a.u.)')




%% 

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fit_wtmkob.pdf
saveas(gca,'fit_wtmkob.fig')
close all; !open fit_wtmkob.pdf
save ./simData/fit_wtmkob.mat 