%% figs1. The figure for same and different processing rate 
set(0,'DefaultAxesColorOrder',[54 100 139; 135 206 255;79 148 205]/255)
% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData = csvread('../expdata/nascent.csv',1,0);
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

plot_flag = 1; 
pars = getParams(); % mko parameters

N = 10;%20-1;
M = 19;
%pr_rate = [linspace(0.1,1,M) linspace(2,10,M-1);];
pr_rate = linspace(1,10,N) ;%linspace(2,10,M-1);];
km_rate = linspace(1,10,M);

score = zeros(N,M); 
for i =1:N % different pr vs same pr. 
    for j = 1:M
    input_pars(1) = pr_rate(i);
    input_pars(2) = km_rate(j);
    residues = calScore_mko(input_pars,nfkb_exp,expData,0);
    score(i,j) = sum(residues.^2);
    end
    disp(i)
end


%csvwrite('./simData/fig2s_wt2.csv',score)
score = sqrt(score/numel(residues));

minimal_score = min(score(:));
figure('units','inch','position',[12 6 8 3])
if plot_flag
    subplot 121
    imagesc(km_rate,pr_rate,score)
    xlabel('K_{mtr}');ylabel('k_{Pr} fold reduction in mko');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 

j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 
subplot 122
calScore_mko([pr_rate(i),km_rate(j)],nfkb_exp,expData,1);

xlim([0 121]);ylim([0 2])
text(30,1.8,strcat('RMSD=',num2str(minimal_score)))
text(30,1.6,strcat('k_{pr}fold_{mko}=',num2str(pr_rate(i)),';Km_{tr}=',num2str(km_rate(i))))
xlabel('Time (min)');ylabel('nascent (a.u.)')

%% 

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fig2s_mko_n2.pdf
saveas(gca,'fig2s_mko_n2.fig')
close all; 
!open fig2s_mko_n2.pdf


save ./simData/fig2s_mko_n2.mat 