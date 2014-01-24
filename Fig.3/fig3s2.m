% path & exp data 
addpath('../src/')
nascent_exp = csvread('../expdata/nascent.csv',1,0);
expData = csvread('../expdata/mRNA.csv',1,0);
nascent_exp(:,[3,5,7])=[]; % delete the std cols 

% plot options 
plot_flag = 1; 

N = 10;
kdegs = linspace(0.01,0.1,N);
kps = linspace(0.1,1,N)
score = zeros(N,N);

for i = 1:N
    for j = 1:N
        
        input_pars(1) = kps(i);
        input_pars(2) = kdegs(j);

        % calculate score 
        
        residues = calScoreCustom(input_pars,nascent_exp,expData);
        score(i,j) = sum(residues.^2)/(numel(residues)-2-1);        
        disp(((i-1)*N+j)/(N*N))
    end
end

save ./simData/fig3s2.mat 
csvwrite('./simData/fig3s2.csv',score)

% write the data into the file 
%csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])
figure
if plot_flag
    imagesc(kdegs,kps,log10(score))%,[-1 2])%,[log10(min(score(:))),2])
    colorbar
    xlabel('k_{deg}^{tko}');ylabel('k_{Pr}^{wt}');
    h = colorbar;
    ylabel(h,'log10(\chi^2)')
end 
