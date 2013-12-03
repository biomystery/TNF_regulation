function dydt = ode(t,x,options,nfkb_input,k)
    nfkb = interp1(nfkb_input(:,1),nfkb_input(:,2),t);
    dydt= 27/6*nfkb^3/(nfkb^3+.5^3) -k*x ;
end
