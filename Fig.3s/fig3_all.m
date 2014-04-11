set(0,'DefaultAxesColorOrder',[54 100 139; 135 206 255;79 148 205]/255)

% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 
expData = csvread('../expdata/mRNA.csv',1,0);
expData = expData(1:4,:); 
expData(1,3) = 0.01;
expData(:,2:end) = expData(:,2:end)/max(expData(:,2));  % normalized expData

% plot options 
plot_flag = 1; 

N = 19;
M = 19;
pr_mko = linspace(1,10,N);
pr_tko = linspace(1,10,N);

score = zeros(N,M); 

for i = 1:N
    for j = 1:M
    input_pars(1) = pr_mko(i);
    input_pars(2) = pr_tko(j); 
    
    residues = calScore(input_pars,nascent_exp,expData,0);
    score(i,j) = sum(residues.^2);
    end 
    disp(i)
end

score = sqrt(score/numel(residues))


% write the data into the file 
%csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])
figure('units','inch','position',[12 6 8 3])
if plot_flag
    subplot 121
    imagesc(pr_tko,pr_mko,score)
    xlabel('k_{Pr} fold reduction in tko');ylabel('k_{pr} fold_{mko}');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 
minimal_score = min(score(:));
j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 

subplot 122
calScore([pr_mko(i) pr_tko(j)],nascent_exp,expData,1);
xlim([0 120]);ylim([0 1.4])
xlabel('Time (min)');ylabel('mRNA (a.u.)')

text(50,1.3,strcat('RMSD=',num2str(minimal_score)))
text(50,1.1,strcat('k_{pr}fold_{mko}=',num2str(pr_mko(i)),',k_{pr}fold_{tko}=',num2str(pr_tko(j))))

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fit_all_2p.pdf
saveas(gca,'fit_all_2p.fig')
close all;!open fit_all_2p.pdf
save ./simData/fig3_all_2p.mat 