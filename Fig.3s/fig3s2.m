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

N = 21;
pr_fold_mko = linspace(1,11,N);
pr_fold_tko = linspace(0.1,3,N+9);
%pr_fold_tko = linspace(1,11,N);
score = zeros(N,N);

for i = 1:N
    for j = 1:N+9
        
        input_pars(1) = pr_fold_mko(i);
        input_pars(2) = pr_fold_tko(j);

        % calculate score 
        
        residues = calScore(input_pars,nascent_exp,expData,0);
        score(i,j) = sum(residues.^2);
        disp(((i-1)*N+j)/(N*N))
    end
end

save ./simData/fig3s2.mat 
csvwrite('./simData/fig3s2.csv',score)

% write the data into the file 
%csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])
figure
if plot_flag
    imagesc(pr_fold_tko,pr_fold_mko,(score),[min(score(:)),min(score(:))+0.5])
    colorbar
    xlabel('k_{Pr}^{tko}');ylabel('k_{Pr}^{mko}');
    h = colorbar;
    ylabel(h,'RMSD')
end 

% find minimal 
minimal_score = min(score(:));
j = find(min(score) == minimal_score);
i = find(score(:,j) == minimal_score); 
calScore([pr_fold_mko(i),pr_fold_tko(j)],nascent_exp,expData,1);


%% run R 

% run R
!R CMD BATCH fig3s2.r
!rm *.Ro* 
