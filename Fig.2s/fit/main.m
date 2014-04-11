% read data from R 
addpath('../../src/')

p = parsFinal;%log10([0.5 0.6 1.5 1.5 3]);
%p = log10([1 1 1 1 1]); 
for i = 1:1
    %p = parsFinal;% + log10(1+(rand(1,5)-0.5)*.4);
    ub = log10([ 1 1 100 100 6]); % V_tr, Km_tr, k_pr, pr_fold
    lb = log10([ 0.05 0.01 0.01 0.01 0.5]);
    ub = p + log10(5) ; lb =p - log10(5); 
    options =optimset('PlotFcns',{@optimplotx,@optimplotfunccount,@optimplotresnorm,@optimplotfirstorderopt},...
                      'TolFun',1e-10,'TolX',1e-10); 

    [parsFinal,resnorm,residual,~,~,~,jacobian] = lsqnonlin(@ ...
                                                      calScore,p,lb, ...
                                                      ub,options);
    chisq = resnorm/(7-5-1); 
    save fitCust.mat

end 