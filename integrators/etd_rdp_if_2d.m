function u = etd_rdp_if_2d(u,A,g,nsteps,tau)

t = 0;
d = 2;
r(1) = 1/3;
r(2) = 1/4;
r(3) = 1;
sz_u = size(u{1});
pn = prod(sz_u);
Id = speye(pn);
for ii= 1:d
  for jj = 1:3
    [L{1}{jj}{ii},U{1}{jj}{ii}]=lu(Id - (r(jj)*tau)*A{1}{ii});
    [L{2}{jj}{ii},U{2}{jj}{ii}]=lu(Id - (r(jj)*tau)*A{2}{ii});
  end
end

u{1} = u{1}(:);
u{2} = u{2}(:);
for jj = 1:nsteps
    gn = g(u{1},u{2});
    gnu = gn(1:pn);
    gnv = gn(pn+1:2*pn);
    % For u
    p1 = U{1}{3}{1}\(L{1}{3}{1}\(u{1} + tau*gnu));
    u_star = U{1}{3}{2}\(L{1}{3}{2}\p1);
    % For v
    p1 = U{2}{3}{1}\(L{2}{3}{1}\(u{2} + tau*gnv));
    v_star = U{2}{3}{2}\(L{2}{3}{2}\p1);

    g_star = g(u_star,v_star);
    g_staru = g_star(1:pn);
    g_starv = g_star(pn+1:2*pn);

    % For u
    b1 = U{1}{1}{1}\(L{1}{1}{1}\gnu);
    b2 = U{1}{2}{1}\(L{1}{2}{1}\gnu);
    c2u = 9*b1-8*b2;
    % For v
    b1 = U{2}{1}{1}\(L{2}{1}{1}\gnv);
    b2 = U{2}{2}{1}\(L{2}{2}{1}\gnv);
    c2v = 9*b1-8*b2;

    % For u
    a1 = U{1}{1}{1}\(L{1}{1}{1}\u{1});
    a2 = U{1}{2}{1}\(L{1}{2}{1}\u{1});
    c1u = 9*a1-8*a2;
    s1u = U{1}{1}{2}\(L{1}{1}{2}\(9*c1u+2*tau*c2u+tau*g_staru));
    s2u = U{1}{2}{2}\(L{1}{2}{2}\(8*c1u+3*tau/2*c2u+tau/2*g_staru));
    u{1} = s1u-s2u;
    % For v
    a1 = U{2}{1}{1}\(L{2}{1}{1}\u{2});
    a2 = U{2}{2}{1}\(L{2}{2}{1}\u{2});
    c1v = 9*a1-8*a2;
    s1v = U{2}{1}{2}\(L{2}{1}{2}\(9*c1v+2*tau*c2v+tau*g_starv));
    s2v = U{2}{2}{2}\(L{2}{2}{2}\(8*c1v+3*tau/2*c2v+tau/2*g_starv));
    u{2} = s1v-s2v;

    t = t + tau;
end
u{1} = reshape(u{1},sz_u);
u{2} = reshape(u{2},sz_u);
end
