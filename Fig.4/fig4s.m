addpath('../src/')

% read data from R 
mRNA_exp = csvread('../expdata/mRNA.csv',1,0);

%% pro TNF
N = 10;
M = 10;
tl_fold = linspace(1,10,N);
sec_fold = linspace(1,10,M);
mRNA_back = mRNA_all; 
rmsd_pro = zeros(N,M);
nrmsd_pro = zeros(N,M);
%
for i = 1:N
    for j = 1:M
        
        input_pars(1) = tl_fold(i);
        input_pars(2) = sec_fold(j);
        
        % calculate residues
        residues = calScoreCustom(input_pars,mRNA_exp);
        score(i,j) = sum(residues.^2)/(numel(residues)-2-1);        
        disp(((i-1)*N+j)/(N*N))

    end
end
save ./simData/fig4s.mat 

%%
figure
if plot_flag
    imagesc(sec_fold,tl_fold,log10(score))
    colorbar
    xlabel('ksec_{fold}');ylabel('ktl_{fold}');
    h = colorbar;
    ylabel(h,'log10(\chi^2)')
end 
