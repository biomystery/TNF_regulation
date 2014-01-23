% read data from R 
addpath('../src/')

nfkb_exp = csvread('../expdata/nfkb.csv',1,0);
expData =  csvread('../expdata/nascent.csv',1,0);
plot_flag = 1; 
pars = getParams(); % wt parameters


N = 31;
kms = linspace(0.3,0.6,N);
kps = linspace(0.04,1.24,N)

score = zeros(N,N);


for i =1:N
    for j = 1:N
        input_pars(1) = kms(i);
        input_pars(2) = kps(j);

        % cal score
        residues = calScoreCustom(input_pars,nfkb_exp,expData);
        
        score(i,j) = sum(residues.^2)/(numel(residues)-2-1);
        
        disp(((i-1)*N+j)/(N*N))
    end
end

save ./simData/fig2s2.mat 

% plot 
if plot_flag
    imagesc(kms,kps,log10(score));%,[0.6 3])
    colorbar
    xlabel('k_{Pr}^{wt}');ylabel('k_{dTr}');
end 

% save data 
csvwrite('./simData/nrmsd2.csv',score)

