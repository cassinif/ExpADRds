function solver_matlab(tstar,Afun,U,g,method,options)
odefun = @(t,u) Afun(u)+g(t,u);
feval(method,odefun,[0,tstar],U,options);
