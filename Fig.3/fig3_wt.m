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

N = 10;
pr_rate = linspace(0.1,1,N);

for i = 1:N
    input_pars(1) = pr_rate(i);
    residues = calScore_wt(input_pars,nascent_exp,expData,0);
    score(i) = sum(residues.^2);
    disp(i)
end

save ./simData/fig3_wt.mat 
csvwrite('./simData/fig3_wt.csv',score)

% write the data into the file 
%csvwrite(filenames{i},[t nascent_wt nascent_mko nascent_tko])
figure
if plot_flag
    plot(pr_rate,score)
    xlabel('k_{Pr}^{wt}');ylabel('RMSD');
end 

% find minimal 
i = find(score == min(score)) ; 

calScore_wt(pr_rate(i),nascent_exp,expData,1);


%% run R 

% run R
%!R CMD BATCH fig3s2.r
%!rm *.Ro* 
