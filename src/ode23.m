function dydt = ode23(t,x,options,nfkb_input,k,kdeg)
    dydt = zeros(2,1);
    kdeg_real = interp1([0 30 60 120],[0.07 kdeg 0.07 0.14],t);    
    nfkb = interp1(nfkb_input(:,1),nfkb_input(:,2),t);
    dydt(1)= 27/6*nfkb^3/(nfkb^3+.5^3) -k*x(1) ;
    dydt(2)= k*x(1) - kdeg_real*x(2) ;    
end
