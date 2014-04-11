function dydt = ode3s(t,x,options,nascent,pars)
    dydt= pars('k_pr')*interp1(nascent(:,1),nascent(:,2),t) - ...
          pars('kdeg_m')*x ;

end
