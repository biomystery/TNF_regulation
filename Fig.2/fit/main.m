% read data from R 
addpath('../../src/')

%p = log10([0.5 0.6 1.5 1.5 3]);

%best_p = parsFinal;
%best_score = resnorm; 

for i = 1:10
           p = parsFinal;% + log10(1+(rand(1,5)-0.5)*.4);
                         %ub = log10([ 1 1 100 100 6]); % V_tr, Km_tr, k_pr, pr_fold
                         %lb = log10([ 0.05 0.01 0.01 0.01 0.5]);
                  ub = p + log10(5) ; lb =p - log10(5); 
    options =optimset('PlotFcns',{@optimplotx,@optimplotfunccount,@optimplotresnorm,@optimplotfirstorderopt},...
                      'TolFun',1e-10,'TolX',1e-10); 

    [parsFinal,resnorm,residual,~,~,~,jacobian] = lsqnonlin(@ ...
                                                      calScore,p,lb, ...
                                                      ub,options);
% $$$     if resnorm < best_score
% $$$         best_score = resnorm;
% $$$         best_p = parsFinal;
% $$$     else
% $$$         parsFinal = best_p; 
% $$$     end
    
    save fit.mat

end 