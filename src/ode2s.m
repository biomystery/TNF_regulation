function dydt = ode2s(t,x,options,nfkb_input,pars)
    nfkb = interp1(nfkb_input(:,1),nfkb_input(:,2),t);
    dydt= pars('V_tr') * nfkb - pars('k_pr')*x ;
end
