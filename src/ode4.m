function dydt = ode4(t,x,options,mRNA,pars)
    dydt= pars('k_tl')*interp1(mRNA(:,1),mRNA(:,2),t) -pars('kdeg_p')*x;
end
