% read data from R 
addpath('../src/')
nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData =  csvread('../expdata/nascent.csv',1,0);
plot_flag = 1; 
pars = getParams(); % wt parameters


N = 10;
M = 21;
%pr_fold = linspace(0.1,3,N);
kd = linspace(0.1,1,N);
kpr = linspace(-2,0,M);
kpr = 10.^kpr; 

score = zeros(N,M);


for i =1:N
    for j = 1:M
        input_pars(1) = kd(i);
        input_pars(2) = kpr(j);

        % cal score
        residues = calScoreCustom(input_pars,nfkb_exp,expData);
        
        score(i,j) = sqrt(sum(residues.^2)/numel(residues));
        
        disp(((i-1)*N+j)/(N*M))
    end
end

save ./simData/fig2s.mat 

% plot 
if plot_flag
    imagesc(kd,kpr,score)
end 

% save data 
csvwrite('./simData/nrmsd.csv',score)

% run R
!R CMD BATCH fig2s.R
!rm *.Ro* 
!rm *.Rh*
!rm *.RD*