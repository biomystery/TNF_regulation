function dydt = ode2_4a(t,x,options,nfkb_input,pars)
    dydt = zeros(3,1);
    kdeg_real = interp1([0 30 60 120],[0.07 pars('kdeg_m') 0.07 0.07],t);    
    nfkb = interp1(nfkb_input(:,1),nfkb_input(:,2),t);
    
    dydt(1)=  pars('V_tr') * nfkb^pars('n') /(nfkb^pars('n') + ...
                                          pars('Km_tr')^pars('n')) - pars('k_pr')*x(1) ;
    dydt(2)= pars('k_pr') * x(1) - kdeg_real*x(2) ;    
    dydt(3)= pars('k_tl') * x(2) - pars('kdeg_p')*x(3) -pars('k_sec')*x(3) ;    
end
