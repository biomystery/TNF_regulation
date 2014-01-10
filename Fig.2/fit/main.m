% read data from R 
addpath('../../src/')

p = log10([11.25 0.5 0.6 1.5 1.5]);

for i = 1:10
    p = parsFinal;
    %ub = log10([100 1 1 100 100]); % V_tr, Km_tr, k_pr, pr_fold
    %lb = log10([0.01 0.05 0.01 0.01 0.01]);
    ub = p + log10(5) ; lb =p - log10(5); 
    options =optimset('PlotFcns',{@optimplotx,@optimplotfunccount,@optimplotresnorm,@optimplotfirstorderopt},...
                      'TolFun',1e-10,'TolX',1e-10); 

    [parsFinal,resnorm,residual,~,~,~,jacobian] = lsqnonlin(@calScore,p,lb,ub,options);

end 