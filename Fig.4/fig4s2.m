addpath('../src/')

% read data from R 
mRNA_exp = csvread('../expdata/mRNA.csv',1,0);

%% pro TNF
N = 10;
M = 10;
%tl_fold = linspace(0.5,10,N);
%sec_fold = linspace(0.5,10,M);

k_degPs = linspace(0.01,0.1,N);
k_secs = linspace(0.1,1,M);

plot_flag = 1; 

%
for i = 1:N
    for j = 1:M
        

        input_pars(1) = k_degPs(i);
        input_pars(2) = k_secs(j);        
        
        % calculate residues
        residues = calScoreCustom(input_pars,mRNA_exp);
        score(i,j) = sum(residues.^2)/(numel(residues)-2-1);        
        disp(((i-1)*N+j)/(N*N))

    end
end
save ./simData/fig4s2.mat 

%%
figure
if plot_flag
    imagesc(k_secs,k_degPs,log10(score))
    colorbar
    xlabel('k_{sec}');ylabel('k_{degP}');
    h = colorbar;
    ylabel(h,'log10(\chi^2)')
end 

