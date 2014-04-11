% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
expData = csvread('../expdata/mRNA.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 

% plot options 
plot_flag = 1; 

N = 21;
pr_fold_mko = linspace(1,11,N);
%pr_fold_tko = linspace(1,3,N);
pr_fold_tko = linspace(0.05,1.05,N);
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

save ./simData/fig3s1_2.mat 
csvwrite('./simData/fig3s1_2.csv',score)

% write the data into the file 
%csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])
figure
if plot_flag
    imagesc(pr_fold_tko,pr_fold_mko,log10(score),[log10(min(score(:))),2])
    colorbar
    xlabel('k_{Pr}^{tko}');ylabel('k_{Pr}^{mko}');
    h = colorbar;
    ylabel(h,'log10(\chi^2)')
end 
