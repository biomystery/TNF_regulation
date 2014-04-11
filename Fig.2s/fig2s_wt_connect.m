filenames = {'./simData/fig2s_wt.mat','./simData/fig2s_wt2.mat','./simData/fig2s_wt3.mat','./simData/fig2s_wt4.mat','./simData/fig2s_wt5.mat','./simData/fig2s_wt6.mat'};     

score_t =[];

for i = 1:6
    disp(i)
    load(filenames{i});
    score_t = [score_t score]; 
end

pr_rate = linspace(0.1,1,10);
km_rate = linspace(0.1,10,100);
minimal_score = min(score_t(:));


imagesc(km_rate,pr_rate,score_t,[minimal_score,minimal_score+0.5])
xlabel('k_{M}');ylabel('k_{Pr}');
h = colorbar;
ylabel(h,'RMSD')


