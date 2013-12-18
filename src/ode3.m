function dydt = ode3(t,x,options,nascent,pars)
    kdeg_real = interp1([0 30 60 120],[0.07 pars('kdeg_m') 0.07 0.07],t);
    dydt= pars('k_pr')*interp1(nascent(:,1),nascent(:,2),t) -kdeg_real*x ;
end