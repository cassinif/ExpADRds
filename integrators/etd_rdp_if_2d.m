function u = etd_rdp_if_2d_thom(u,A,g,nsteps,tau)

t = 0;
d = 2;
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
  M{1}{jj}{1}=Id - (r(jj)*tau)*A{1}{1};
  ldiag{1}{jj}{1} = full([0;diag(M{1}{jj}{1},-1)]);
  ddiag{1}{jj}{1} = full(diag(M{1}{jj}{1}));
  udiag{1}{jj}{1} = full([diag(M{1}{jj}{1},1);0]);
  
  M{2}{jj}{1}=Id - (r(jj)*tau)*A{2}{1};
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
u{1} = u{1}(:);
u{2} = u{2}(:);

for jj = 1:nsteps
    gn = g(u{1},u{2});
    gnu = gn(1:pn);
    gnv = gn(pn+1:2*pn);
    % For u
    p1 = thom_ad(ldiag{1}{3}{1},ddiag{1}{3}{1},udiag{1}{3}{1},1,u{1}+tau*gnu);
    u_star = thom_ad(ldiag{1}{3}{2},ddiag{1}{3}{2},udiag{1}{3}{2},sz_u(1),p1);
    % For v
    p1 = thom_ad(ldiag{2}{3}{1},ddiag{2}{3}{1},udiag{2}{3}{1},1,u{2}+tau*gnv);
    v_star = thom_ad(ldiag{2}{3}{2},ddiag{2}{3}{2},udiag{2}{3}{2},sz_u(1),p1);

    g_star = g(u_star,v_star);
    g_staru = g_star(1:pn);
    g_starv = g_star(pn+1:2*pn);

    % For u
    b1 = thom_ad(ldiag{1}{1}{1},ddiag{1}{1}{1},udiag{1}{1}{1},1,gnu);
    b2 = thom_ad(ldiag{1}{2}{1},ddiag{1}{2}{1},udiag{1}{2}{1},1,gnu);
    c2u = 9*b1-8*b2;
    % For v
    b1 = thom_ad(ldiag{2}{1}{1},ddiag{2}{1}{1},udiag{2}{1}{1},1,gnv);
    b2 = thom_ad(ldiag{2}{2}{1},ddiag{2}{2}{1},udiag{2}{2}{1},1,gnv);
    c2v = 9*b1-8*b2;

    % For u
    a1 = thom_ad(ldiag{1}{1}{1},ddiag{1}{1}{1},udiag{1}{1}{1},1,u{1});
    a2 = thom_ad(ldiag{1}{2}{1},ddiag{1}{2}{1},udiag{1}{2}{1},1,u{1});
    c1u = 9*a1-8*a2;
    s1u = thom_ad(ldiag{1}{1}{2},ddiag{1}{1}{2},udiag{1}{1}{2},sz_u(1),9*c1u+2*tau*c2u+tau*g_staru);
    s2u = thom_ad(ldiag{1}{2}{2},ddiag{1}{2}{2},udiag{1}{2}{2},sz_u(1),8*c1u+3*tau/2*c2u+tau/2*g_staru);
    u{1} = s1u-s2u;
    % For v
    a1 = thom_ad(ldiag{2}{1}{1},ddiag{2}{1}{1},udiag{2}{1}{1},1,u{2});
    a2 = thom_ad(ldiag{2}{2}{1},ddiag{2}{2}{1},udiag{2}{2}{1},1,u{2});
    c1v = 9*a1-8*a2;
    s1v = thom_ad(ldiag{2}{1}{2},ddiag{2}{1}{2},udiag{2}{1}{2},sz_u(1),9*c1v+2*tau*c2v+tau*g_starv);
    s2v = thom_ad(ldiag{2}{2}{2},ddiag{2}{2}{2},udiag{2}{2}{2},sz_u(1),8*c1v+3*tau/2*c2v+tau/2*g_starv);
    u{2} = s1v-s2v;

    t = t + tau;
end

u{1} = reshape(u{1},sz_u);
u{2} = reshape(u{2},sz_u);
end
