function dydt = ode2(t,x,options,mRNA,kdeg,ksec)
dydt= interp1(mRNA(:,1),mRNA(:,2),t) -kdeg*x -ksec*x ;
end
