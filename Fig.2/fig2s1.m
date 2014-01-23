% read data from R 
addpath('../src/')

nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData =  csvread('../expdata/nascent.csv',1,0);
plot_flag = 1; 
pars = getParams(); % wt parameters


N = 30;
pr_fold_mko = linspace(3,6,N);
pr_fold_tko = linspace(1,4,N);
score = zeros(N,N);


for i =1:N
    for j = 1:N
        input_pars(1) = pr_fold_mko(i);
        input_pars(2) = pr_fold_tko(j);

        % cal score
        residues = calScoreCustom(input_pars,nfkb_exp,expData);
        
        score(i,j) = sum(residues.^2)/(numel(residues)-2-1);
        
        disp(((i-1)*N+j)/(N*N))
    end
end

save ./simData/fig2s.mat 

% plot 
if plot_flag
    imagesc(pr_fold_tko,pr_fold_mko,score,[0.6 3])
    colorbar
    xlabel('k_{Pr}^{tko}');ylabel('k_{Pr}^{mko}');
end 

% save data 
csvwrite('./simData/nrmsd.csv',score)

% run R
!R CMD BATCH fig2s.R
!rm *.Ro* 
