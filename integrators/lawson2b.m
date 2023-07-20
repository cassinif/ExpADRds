function U = lawson2b(U,A,g,nsteps,tau)
t = 0;
d = length(A{1});
for mu = 1:d
  E{1}{mu} = expm(tau*A{1}{mu});
  E{2}{mu} = expm(tau*A{2}{mu});
end
for jj = 1:nsteps
  gnu = g{1}(t,U{1},U{2});
  gnv = g{2}(t,U{1},U{2});

  U2{1} = tucker(U{1} + tau*gnu,E{1});
  U2{2} = tucker(U{2} + tau*gnv,E{2});
  U{1} = tucker(U{1}+tau/2*gnu,E{1})+tau/2*g{1}(t+tau,U2{1},U2{2});
  U{2} = tucker(U{2}+tau/2*gnv,E{2})+tau/2*g{2}(t+tau,U2{1},U2{2});
  t = t + tau;
end
end
