% read data from R 
addpath('../src/')

p = log10([0.5 3]);

for i = 1:10
    % p = parsFinal;
    lb = log10([ 0.1 0.5]);
    ub = log10([1 6]);
    %ub = p + log10(5) ; lb =p - log10(5); 
    options =optimset('PlotFcns',{@optimplotx,@optimplotfunccount,@optimplotresnorm,@optimplotfirstorderopt},...
                      'TolFun',1e-10,'TolX',1e-10); 

    [parsFinal,resnorm,residual,~,~,~,jacobian] = lsqnonlin(@ ...
                                                      calScoreCustom23,p,lb, ...
                                                      ub,options);
    chisq = resnorm/(numel(residual)-numel(p)-1); 
    save fitCust23.mat

end 
