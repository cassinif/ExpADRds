function U = etd2rkds(U,A,F,g,nsteps,tau)
  t = 0;
  d = length(A{1});
  for mu = 1:d
    [P1{1}{mu},P2{1}{mu}] = phiquad(tau*A{1}{mu},2);
    [P1{2}{mu},P2{2}{mu}] = phiquad(tau*A{2}{mu},2);
  end
  for jj = 1:nsteps
    U2{1} = U{1} + tau*tucker(F{1}(t,U{1},U{2}),P1{1});
    U2{2} = U{2} + tau*tucker(F{2}(t,U{1},U{2}),P1{2});
    Up{1} = U2{1} + (2^(d-1)*tau)*tucker(g{1}(t+tau,U2{1},U2{2})-g{1}(t,U{1},U{2}),P2{1});
    U{2} = U2{2} + (2^(d-1)*tau)*tucker(g{2}(t+tau,U2{1},U2{2})-g{2}(t,U{1},U{2}),P2{2});
    U{1} = Up{1};
    t = t + tau;
  end
end
