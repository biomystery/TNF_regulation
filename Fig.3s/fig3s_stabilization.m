% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
expData = csvread('../expdata/mRNA.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 

% plot options 
plot_flag = 1; 

N = 10;
pr_fold_mko = linspace(.01,.1,N);
%pr_fold_tko = linspace(1,3,N);
pr_fold_tko = linspace(0.01,.1,N);
score = zeros(N,N);

for i = 1:N
    for j = 1:N
        
        input_pars(1) = pr_fold_mko(i);
        input_pars(2) = pr_fold_tko(j);

        % calculate score 
        
        residues = calScoreCustom(input_pars,nascent_exp,expData);
        score(i,j) = sum(residues.^2)/(numel(residues)-2-1);        
        disp(((i-1)*N+j)/(N*N))
    end
end

save ./simData/fig3s_stabilization.mat 


figure
if plot_flag
    imagesc(pr_fold_tko,pr_fold_mko,(score))
    colorbar
    xlabel('kdeg_{m}^{tko}');ylabel('kdeg_{m}^{wt/mko}');
    h = colorbar;
    ylabel(h,'log10(\chi^2)')
    set(gca,'yDir','normal')
end 
