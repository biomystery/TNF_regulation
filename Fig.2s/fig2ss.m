%% figs1. The figure for same and different processing rate 
set(0,'DefaultAxesColorOrder',[54 100 139; 135 206 255;79 148 205]/255)
% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData = csvread('../expdata/nascent.csv',1,0);
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

plot_flag = 1; 
pars = getParams(); % wt parameters
k_pr = pars('k_pr');

N = 10;
M = 10;
pr_rate = linspace(1,10,N);
km_rate = linspace(1.1,2,M);

score = zeros(N,M); 
for i =1:N % different pr vs same pr. 
    for j = 1:M
    input_pars(1) = pr_rate(i);
    input_pars(2) = km_rate(j);
    residues = calScores(input_pars,nfkb_exp,expData,0);
    score(i,j) = sum(residues.^2);
    end
    disp(i)
end

score = sqrt(score/numel(residues));


%csvwrite('./simData/fig2s_wt2.csv',score)


minimal_score = min(score(:));
figure('units','inch','position',[12 6 8 3])
if plot_flag
    subplot 121
    imagesc(km_rate,pr_rate,score)
    xlabel('k_{Mpr}fold_{tko}');ylabel('k_{MPr}fold_{mko}');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 

j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 

subplot 122
calScores([pr_rate(i(1)),km_rate(j(1))],nfkb_exp,expData,1);
text(50,1.8,strcat('RMSD=',num2str(score(i,j))))
text(50,1.6,strcat('k_{mpr}fold_{mko}=',num2str(pr_rate(i)),',k_{mpr}fold_{tko}=',num2str(km_rate(j))))

xlim([0 121]);ylim([0 3])
xlabel('Time (min)');ylabel('nascent (a.u.)')

%% 

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fitss_all_km.pdf
saveas(gca,'fitss_all_km.fig')
close all; !open fitss_all_km.pdf
save ./simData/fitss_all_km2.mat 