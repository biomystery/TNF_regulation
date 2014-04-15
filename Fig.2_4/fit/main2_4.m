clear;
set(0,'DefaultAxesColorOrder',[54 100 139; 135 206 255;79 148 205]/255)
% read data from R 
addpath('../../src/')

parMatirx = [
    0.65 0.001 10 ;% km_tr
    0.4 0.1 1; % k_pr
    4 2 5; % fold_pr_mko 
    1.5 1 2; %fold_pr_tko
    2 1 3; % fold_kmtr_mko 
    1 0.001 100; % tl 
    3.5 1 5; %fold_tl_tko
    5 1 10; % fold_sec_tko 
    1 .001 100; % k_sec
    0.5 0.001 1; % kdeg_p
    ];

pstart = log10(parMatirx(:,1));
lb = log10(parMatirx(:,2));
ub = log10(parMatirx(:,3));


options = saoptimset('PlotFcns',{@saplotbestx,@saplottemperature,@ ...
                    saplotf,@saplotstopping,@saplotbestf});

options = saoptimset(options,'TemperatureFcn',@temperaturefast);
options = saoptimset(options,'ReannealInterval',50);
options = saoptimset(options,'Display','iter','DisplayInterval',400);
options = saoptimset(options,'TolFun',1e-5);

[x,fval,exitFlag,output] = simulannealbnd(@objectFun, pstart,lb,ub, ...
                                          options)

fprintf('The number of iterations was : %d\n', output.iterations);
fprintf('The number of function evaluations was : %d\n', output.funccount);
fprintf('The best function value found was : %g\n', fval);


save optim.mat
