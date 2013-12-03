function dydt = ode(t,x,options,mRNA,kdeg)
dydt= interp1(mRNA(:,1),mRNA(:,2),t) -kdeg*x ;
end
