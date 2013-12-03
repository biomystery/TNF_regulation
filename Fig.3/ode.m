function dydt = ode(t,x,options,nascent,kdeg)
    kdeg_real = interp1([0 30 60 120],[0.07 kdeg 0.07 0.07],t);
    dydt= interp1(nascent(:,1),nascent(:,2),t) -kdeg_real*x ;
end
