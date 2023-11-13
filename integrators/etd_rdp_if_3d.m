function u = etd_rdp_if_3d(u,A,g,nsteps,tau)

t = 0;
d = 3;
r(1) = 1/3;
r(2) = 1/4;
r(3) = 1;

sz_u = size(u{1});
pn = prod(sz_u);
Id = speye(pn);

% We could avoid creating the big sparse matrices, but negligible cost
% when time integrating
% x direction
for jj = 1:3
  M{1}{jj}{1}= Id - (r(jj)*tau)*A{1}{1};
  ldiag{1}{jj}{1} = full([0;diag(M{1}{jj}{1},-1)]);
  ddiag{1}{jj}{1} = full(diag(M{1}{jj}{1}));
  udiag{1}{jj}{1} = full([diag(M{1}{jj}{1},1);0]);


  M{2}{jj}{1}= Id - (r(jj)*tau)*A{2}{1};
  ldiag{2}{jj}{1} = full([0;diag(M{2}{jj}{1},-1)]);
  ddiag{2}{jj}{1} = full(diag(M{2}{jj}{1}));
  udiag{2}{jj}{1} = full([diag(M{2}{jj}{1},1);0]);
end

% y direction
for jj = 1:3
  M{1}{jj}{2}=Id - (r(jj)*tau)*A{1}{2};
  ldiag{1}{jj}{2} = full([zeros(sz_u(1),1);diag(M{1}{jj}{2},-sz_u(1))]);
  ddiag{1}{jj}{2} = full(diag(M{1}{jj}{2}));
  udiag{1}{jj}{2} = full([diag(M{1}{jj}{2},sz_u(1));zeros(sz_u(1),1)]);

  M{2}{jj}{2}=Id - (r(jj)*tau)*A{2}{2};
  ldiag{2}{jj}{2} = full([zeros(sz_u(1),1);diag(M{2}{jj}{2},-sz_u(1))]);
  ddiag{2}{jj}{2} = full(diag(M{2}{jj}{2}));
  udiag{2}{jj}{2} = full([diag(M{2}{jj}{2},sz_u(1));zeros(sz_u(1),1)]);
end

% z direction
for jj = 1:3
  M{1}{jj}{3}=Id - (r(jj)*tau)*A{1}{3};
  ldiag{1}{jj}{3} = full([zeros(sz_u(1)*sz_u(2),1);diag(M{1}{jj}{3},-sz_u(1)*sz_u(2))]);
  ddiag{1}{jj}{3} = full(diag(M{1}{jj}{3}));
  udiag{1}{jj}{3} = full([diag(M{1}{jj}{3},sz_u(1)*sz_u(2));zeros(sz_u(1)*sz_u(2),1)]);

  M{2}{jj}{3}=Id - (r(jj)*tau)*A{2}{3};
  ldiag{2}{jj}{3} = full([zeros(sz_u(1)*sz_u(2),1);diag(M{2}{jj}{3},-sz_u(1)*sz_u(2))]);
  ddiag{2}{jj}{3} = full(diag(M{2}{jj}{3}));
  udiag{2}{jj}{3} = full([diag(M{2}{jj}{3},sz_u(1)*sz_u(2));zeros(sz_u(1)*sz_u(2),1)]);
end

u{1} = u{1}(:);
u{2} = u{2}(:);

for jj = 1:nsteps
  gn = g(u{1},u{2});
  gnu = gn(1:pn);
  gnv = gn(pn+1:2*pn);

  % For u
  p1 = thom_ad(ldiag{1}{3}{1},ddiag{1}{3}{1},udiag{1}{3}{1},1,u{1}+tau*gnu);
  p2 = thom_ad(ldiag{1}{3}{2},ddiag{1}{3}{2},udiag{1}{3}{2},sz_u(1),p1);
  u_star = thom_ad(ldiag{1}{3}{3},ddiag{1}{3}{3},udiag{1}{3}{3},sz_u(1)*sz_u(2),p2);
  % For v
  p1 = thom_ad(ldiag{2}{3}{1},ddiag{2}{3}{1},udiag{2}{3}{1},1,u{2}+tau*gnv);
  p2 = thom_ad(ldiag{2}{3}{2},ddiag{2}{3}{2},udiag{2}{3}{2},sz_u(1),p1);
  v_star = thom_ad(ldiag{2}{3}{3},ddiag{2}{3}{3},udiag{2}{3}{3},sz_u(1)*sz_u(2),p2);

  g_star = g(u_star,v_star);
  g_staru = g_star(1:pn);
  g_starv = g_star(pn+1:2*pn);

  % For u
  b1 = thom_ad(ldiag{1}{1}{1},ddiag{1}{1}{1},udiag{1}{1}{1},1,gnu);
  b2 = thom_ad(ldiag{1}{2}{1},ddiag{1}{2}{1},udiag{1}{2}{1},1,gnu);
  c2 = 9*b1-8*b2;
  b3 = thom_ad(ldiag{1}{1}{2},ddiag{1}{1}{2},udiag{1}{1}{2},sz_u(1),c2);
  b4 = thom_ad(ldiag{1}{2}{2},ddiag{1}{2}{2},udiag{1}{2}{2},sz_u(1),c2);
  c4u = 9*b3-8*b4;
  % For v
  b1 = thom_ad(ldiag{2}{1}{1},ddiag{2}{1}{1},udiag{2}{1}{1},1,gnv);
  b2 = thom_ad(ldiag{2}{2}{1},ddiag{2}{2}{1},udiag{2}{2}{1},1,gnv);
  c2 = 9*b1-8*b2;
  b3 = thom_ad(ldiag{2}{1}{2},ddiag{2}{1}{2},udiag{2}{1}{2},sz_u(1),c2);
  b4 = thom_ad(ldiag{2}{2}{2},ddiag{2}{2}{2},udiag{2}{2}{2},sz_u(1),c2);
  c4v = 9*b3-8*b4;

  % For u
  a1 = thom_ad(ldiag{1}{1}{1},ddiag{1}{1}{1},udiag{1}{1}{1},1,u{1});
  a2 = thom_ad(ldiag{1}{2}{1},ddiag{1}{2}{1},udiag{1}{2}{1},1,u{1});
  c1 = 9*a1-8*a2;
  a3 = thom_ad(ldiag{1}{1}{2},ddiag{1}{1}{2},udiag{1}{1}{2},sz_u(1),c1);
  a4 = thom_ad(ldiag{1}{2}{2},ddiag{1}{2}{2},udiag{1}{2}{2},sz_u(1),c1);
  c3u = 9*a3-8*a4;
  s1u = thom_ad(ldiag{1}{1}{3},ddiag{1}{1}{3},udiag{1}{1}{3},sz_u(1)*sz_u(2),9*c3u+2*tau*c4u+tau*g_staru);
  s2u = thom_ad(ldiag{1}{2}{3},ddiag{1}{2}{3},udiag{1}{2}{3},sz_u(1)*sz_u(2),8*c3u+3*tau/2*c4u+tau/2*g_staru);
  u{1} = s1u-s2u;
  % For v
  a1 = thom_ad(ldiag{2}{1}{1},ddiag{2}{1}{1},udiag{2}{1}{1},1,u{2});
  a2 = thom_ad(ldiag{2}{2}{1},ddiag{2}{2}{1},udiag{2}{2}{1},1,u{2});
  c1 = 9*a1-8*a2;
  a3 = thom_ad(ldiag{2}{1}{2},ddiag{2}{1}{2},udiag{2}{1}{2},sz_u(1),c1);
  a4 = thom_ad(ldiag{2}{2}{2},ddiag{2}{2}{2},udiag{2}{2}{2},sz_u(1),c1);
  c3v = 9*a3-8*a4;
  s1v = thom_ad(ldiag{2}{1}{3},ddiag{2}{1}{3},udiag{2}{1}{3},sz_u(1)*sz_u(2),9*c3v+2*tau*c4v+tau*g_starv);
  s2v = thom_ad(ldiag{2}{2}{3},ddiag{2}{2}{3},udiag{2}{2}{3},sz_u(1)*sz_u(2),8*c3v+3*tau/2*c4v+tau/2*g_starv);
  u{2} = s1v-s2v;

  t = t + tau;
end


u{1} = reshape(u{1},sz_u);
u{2} = reshape(u{2},sz_u);
end
