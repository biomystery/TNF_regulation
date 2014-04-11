%% figs1. The figure for same and different processing rate 

% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData = csvread('../expdata/nascent.csv',1,0);
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

plot_flag = 1; 
pars = getParams(); % mko parameters
k_pr = pars('k_pr');

N = 20-1;
M = 10;
pr_rate = [linspace(0.1,1,M) linspace(2,10,M-1);];
km_rate = linspace(1,10,M);

score = zeros(N,M); 
for i =1:N % different pr vs same pr. 
    for j = 1:M
    input_pars(1) = pr_rate(i);
    input_pars(2) = km_rate(j);
    residues = calScore_mkos(input_pars,nfkb_exp,expData,0);
    score(i,j) = sum(residues.^2);
    end
    disp(i)
end

save ./simData/fig2s_mko.mat 
%csvwrite('./simData/fig2s_wt2.csv',score)


minimal_score = min(score(:));
figure
if plot_flag
    imagesc(km_rate,log10(pr_rate),score)
    xlabel('V_{tr}');ylabel('k_{Pr}');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 

j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 
calScore_mkos([pr_rate(i(1)),km_rate(j(2))],nfkb_exp,expData,1);

%% run R 

subplot 122
calScore_mkos([pr_rate(i),km_rate(j)],nfkb_exp,expData,1);

xlim([0 121]);ylim([0 2])
xlabel('Time (min)');ylabel('nascent (a.u.)')

%% 

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fitss_tko_diffpr.pdf
saveas(gca,'fitss_tko_diffpr.fig')