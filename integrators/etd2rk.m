function U = etd2rk_phiks(U,A,F,g,nsteps,tau,tol)
t = 0;
for jj = 1:nsteps
  W{1} = phiks(tau,A{1},F{1}(t,U{1},U{2}),1,tol);
  W{2} = phiks(tau,A{2},F{2}(t,U{1},U{2}),1,tol);
  U2{1} = U{1} + tau*W{1}{2};
  U2{2} = U{2} + tau*W{2}{2};

  W{1} = phiks(tau,A{1},g{1}(t+tau,U2{1},U2{2})-g{1}(t,U{1},U{2}),2,tol);
  W{2} = phiks(tau,A{2},g{2}(t+tau,U2{1},U2{2})-g{2}(t,U{1},U{2}),2,tol);
  U{1} = U2{1} + tau*W{1}{3};
  U{2} = U2{2} + tau*W{2}{3};
  t = t + tau;
end
end
