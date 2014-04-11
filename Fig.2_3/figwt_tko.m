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

N = 10;
M = 40;
pr_rate = linspace(.1,1,N);
km_rate = linspace(1.1,5,M);

score = zeros(N,M); 
for i =1:N % different pr vs same pr. 
    for j = 1:M
    input_pars(1) = pr_rate(i);
    input_pars(2) = km_rate(j);
    residues = calScores_wttko(input_pars,nfkb_exp,expData1,0);
    residues = [residues;calScore_wttko(input_pars,nascent_exp,expData2,0)];
    score(i,j) = sum(residues.^2);
    end
    disp(i)
end

score = sqrt(score/numel(residues));


minimal_score = min(score(:));
figure('units','inch','position',[12 6 12 3])
if plot_flag
    subplot 131
    imagesc(km_rate,pr_rate,score,[minimal_score,minimal_score+.2])
    xlabel('k_{pr}fold_{tko}');ylabel('K_{mtr}');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 

j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 

subplot 132
calScores_wttko([pr_rate(i(1)),km_rate(j(1))],nfkb_exp,expData1,1);
text(50,1.8,strcat('RMSD=',num2str(score(i,j))))
text(50,1.6,strcat('K_{mtr}=',num2str(pr_rate(i)),',k_{pr}fold_{tko}=',num2str(km_rate(j))))

xlim([0 121]);ylim([0 2])
xlabel('Time (min)');ylabel('nascent (a.u.)')

subplot 133
calScore_wttko([pr_rate(i) km_rate(j)],nascent_exp,expData2,1);
xlim([0 120]);ylim([0 1.4])
xlabel('Time (min)');ylabel('mRNA (a.u.)')

text(50,1.3,strcat('RMSD=',num2str(minimal_score)))
text(50,1.1,strcat('k_{mtr}=',num2str(pr_rate(i)),',k_{pr}fold_{tko}=',num2str(km_rate(j))))


%% 

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fit_wttko2.pdf
saveas(gca,'fit_wttko2.fig')
close all; !open fit_wttko2.pdf
save ./simData/fit_wttko2.mat 