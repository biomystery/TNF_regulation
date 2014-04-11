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
pr_rate = [linspace(0.01,.1,10) linspace(0.2,1,9)];
km_rate = linspace(1,10,M);

score = zeros(N,M); 

for i =1:N % different pr vs same pr. 
    for j = 1:M
    input_pars(1) = pr_rate(i);
    input_pars(2) = km_rate(j);
    residues = calScore_wt(input_pars,nfkb_exp,expData,0);
    score(i,j) = sum(residues.^2);
    end
    disp(i)
end

score = sqrt(score/numel(residues)); 


minimal_score = min(score(:));
figure('units','inch','position',[12 6 8 3])
subplot 121
if plot_flag
    imagesc(km_rate,log10(pr_rate),score)
    xlabel('k_{M}');ylabel('log10(k_{Pr})');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 
subplot 122
j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 
calScore_wt([pr_rate(i),km_rate(j)],nfkb_exp,expData,1);
xlim([0 121]);ylim([0 2])
xlabel('Time (min)');ylabel('nascent (a.u.)')
text(50,1.5,strcat('RMSD=',num2str(minimal_score)))
text(50,1.2,strcat('k_{pr}=',num2str(pr_rate(i)),',Km_{tr}=',num2str(km_rate(j))))



set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fit_wt_n2.pdf
saveas(gca,'fit_wt_n2.fig')
close all; 
!open fit_wt_n2.pdf

save ./simData/fig2s_wt_n2.mat 