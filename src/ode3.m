function dydt = ode3(t,x,options,nascent,pars)
    kdeg_real = interp1([0 30 60 120],[pars('kdeg_m_basal') pars('kdeg_m') ...
                        pars('kdeg_m_basal') pars('kdeg_m_basal')],t);
    dydt= pars('k_pr')*interp1(nascent(:,1),nascent(:,2),t) - ...
          kdeg_real*x ;

end
