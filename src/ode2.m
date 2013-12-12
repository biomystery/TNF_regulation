function dydt = ode2(t,x,options,nfkb_input,pars)
    nfkb = interp1(nfkb_input(:,1),nfkb_input(:,2),t);
    dydt= pars('V_tr') * nfkb^pars('n') /(nfkb^pars('n') + ...
                                          pars('Km_tr')^pars('n')) - pars('k_pr')*x ;
end
