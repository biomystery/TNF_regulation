%% figs1. The figure for same and different processing rate 

% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData = csvread('../expdata/nascent.csv',1,0);
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

plot_flag = 1; 
pars = getParams(); % wt parameters
k_pr = pars('k_pr');

N = 19;
M = 10;
pr_rate = linspace(1,10,N);
km_rate = linspace(1,10,M);

score = zeros(N,M); 
for i =1:N % different pr vs same pr. 
    for j = 1:M
    input_pars(1) = pr_rate(i);
    input_pars(2) = km_rate(j);
    residues = calScore_tko(input_pars,nfkb_exp,expData,0);
    score(i,j) = sum(residues.^2);
    end
    disp(i)
end


score = sqrt(score/numel(residues));
minimal_score = min(score(:));
figure('units','inch','position',[12 6 8 3])
subplot 121
if plot_flag
    imagesc(km_rate,pr_rate,score)%,[minimal_score,.2])
    xlabel('k_{M}');ylabel('k_{Pr} fold reduction in tko');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 
subplot 122
j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 
calScore_tko([pr_rate(i(1)),km_rate(j(1))],nfkb_exp,expData,1);
xlim([0 121]);ylim([0 2])
xlabel('Time (min)');ylabel('nascent (a.u.)')
text(50,1.5,strcat('RMSD=',num2str(minimal_score)))
text(50,1.2,strcat('k_{pr}fold_{tko}=',num2str(pr_rate(i)),',Km_{tr}=',num2str(km_rate(j))))
%% run R 

% run R
%!R CMD BATCH fig3s2.r
%!rm *.Ro* 

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fit_same_pr.pdf
saveas(gca,'fit_same_pr.fig')
close all; 
!open fit_same_pr.pdf

save ./simData/fig2s_tko_n2.mat 