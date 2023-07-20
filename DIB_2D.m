clear all
close all

addpath('integrators')
addpath('matfiles')
addpath('external/phisplit')
addpath('external/phisplit/extern/KronPACK/src')
addpath('external/phisplit/extern/phiks')

d = 2;

n = 200*ones(1,d);
a = 0*ones(1,d);
b = 20*ones(1,d);
T = 2.5;

deltau = 1;
rho = 25/4;
a1u = 10;
a2u = 1;
a3u = 66;
a4u = 0.5;
deltav = 20;
a1v = 3;
a2v = 2.5;
a3v = 0.2;
a4v = a1v*(1-a4u)*(1-a3v+a3v*a4u)/a4u/(1+a3v*a4u);
a5v = 1.5;

method = 'ETD2RKds'; % ETD2RKds, ETD2RK, Lawson2b, ETD-RDP-IF, DIRK23, RK32
compute_err = true; % if true, measure error against precomputed reference solution
                     % else plot u component at time T
tol_matlab = 5e-5; % needed if using matlab ODE suite
tol_phiks = 1e-3; % needed if using ETD2RK
nsteps = 1250; % needed if not using matlab ODE suite
tau = T/nsteps;

for mu = 1:d
  x{mu} = linspace(a(mu),b(mu),n(mu));
  h(mu) = (b(mu)-a(mu))/(n(mu)-1);
  D2{mu} = spdiags(ones(n(mu),1)*([1,-2,1]/(h(mu)^2)),-1:1,n(mu),n(mu));
  D2{mu}(1,1:2) = [-2,2]/(h(mu)^2);
  D2{mu}(n(mu),(n(mu)-1):n(mu)) = [2,-2]/(h(mu)^2);
  A_sp{1}{mu} = deltau*D2{mu};
  A_sp{2}{mu} = deltav*D2{mu};
  A{1}{mu} = full(A_sp{1}{mu});
  A{2}{mu} = full(A_sp{2}{mu});
end
[X{1:d}] = ndgrid(x{1:d});

g{1} = @(t,u,v) rho*(a1u*(1-v).*u - a2u*(u.*u).*u - a3u*(v-a4u));
g{2} = @(t,u,v) rho*(a1v*(1+a2v*u).*(1-v).*(1-a3v*(1-v))-a4v*v.*(1+a5v*u).*(1+a3v*v));

dgdu{1}{1} = @(t,u,v) rho*(a1u*(1-v)-(3*a2u)*(u.*u)); %dg1du
dgdu{1}{2} = @(t,u,v) rho*(-a1u*u-a3u); %dg1dv
dgdu{2}{1} = @(t,u,v) rho*((a1v*a2v)*(1-v).*(1-a3v*(1-v))-(a5v*a4v)*v.*(1+a3v*v)); %dg2du
dgdu{2}{2} = @(t,u,v) rho*(a1v*(1+a2v*u).*(-(1-a3v*(1-v))+a3v*(1-v))-a4v*(1+a5v*u).*((1+a3v*v)+a3v*v)); %dg2dv

F{1} = @(t,u,v) kronsumv(u,A{1}) + g{1}(t,u,v);
F{2} = @(t,u,v) kronsumv(v,A{2}) + g{2}(t,u,v);

pn = prod(n);

% For Matlab solver
K{1} = kronsum(A_sp{1});
K{2} = kronsum(A_sp{2});
options.Jacobian = @(t,uvec) [K{1}+spdiags(dgdu{1}{1}(t,uvec(1:pn),uvec(pn+1:2*pn)),0,pn,pn),...
 spdiags(dgdu{1}{2}(t,uvec(1:pn),uvec(pn+1:2*pn)),0,pn,pn);...
 spdiags(dgdu{2}{1}(t,uvec(1:pn),uvec(pn+1:2*pn)),0,pn,pn),...
 K{2}+spdiags(dgdu{2}{2}(t,uvec(1:pn),uvec(pn+1:2*pn)),0,pn,pn)]; % only for DIRK23

options.RelTol = tol_matlab;
options.AbsTol = tol_matlab;
savestr = ['matfiles/',method,'sol'];
options.OutputFcn = @(t,u,flag) myoutfcn(t,u,flag,T,savestr);

Kfun = @(uvec) [K{1}*uvec(1:pn);K{2}*uvec(pn+1:2*pn)];

% For etd_rdp_if
A_otimes{1}{1} = kron(speye(n(2)),A_sp{1}{1});
A_otimes{1}{2} = kron(A_sp{1}{2},speye(n(1)));
A_otimes{2}{1} = kron(speye(n(2)),A_sp{2}{1});
A_otimes{2}{2} = kron(A_sp{2}{2},speye(n(1)));

gvec = @(t,uvec) [g{1}(t,uvec(1:pn),uvec(pn+1:2*pn));g{2}(t,uvec(1:pn),uvec(pn+1:2*pn))];
g_if = @(u,v) gvec(NaN,[u(:);v(:)]);

load('DIB_2D_U0.mat')
u0 = [U0{1}(:);U0{2}(:)];

fprintf('Method: %s\n',method)
switch method
  case 'ETD2RKds'
    tic
    U = etd2rkds(U0,A,F,g,nsteps,tau);
    wctime = toc;
  case 'ETD2RK'
    tic
    U = etd2rk(U0,A,F,g,nsteps,tau,tol_phiks);
    wctime = toc;
  case 'Lawson2b'
    tic
    U = lawson2b(U0,A,g,nsteps,tau);
    wctime = toc;
  case 'ETD-RDP-IF'
    tic
    U = etd_rdp_if_2d(U0,A_otimes,g_if,nsteps,tau);
    wctime = toc;
  case 'DIRK23'
    tic
    solver_matlab(T,Kfun,u0,gvec,'ode23tb',options);
    wctime = toc;
    load(savestr)
    U{1} = reshape(app(1:pn),n);
    U{2} = reshape(app(pn+1:2*pn),n);
  case 'RK32'
    tic
    solver_matlab(T,Kfun,u0,gvec,'ode23',options);
    wctime = toc;
    load(savestr)
    U{1} = reshape(app(1:pn),n);
    U{2} = reshape(app(pn+1:2*pn),n);
  otherwise
    error('Method not known.')
end

if compute_err
  load('DIB_2D_Uref.mat')
  normrifu = norm(Uref{1},'fro');
  normrifv = norm(Uref{2},'fro');
  err = norm([norm(U{1}-Uref{1},'fro')/normrifu,norm(U{2}-Uref{2},'fro')/normrifv]);
  fprintf('Error: %.3e\n',err)
else
  figure;
  surf(X{1},X{2},U{1},'edgecolor','none');
  axis equal
  view(2)
  xlabel('x_1')
  ylabel('x_2')
  colorbar
  drawnow
end

fprintf('Wall-clock time: %.2f s\n',wctime)

rmpath('integrators')
rmpath('matfiles')
rmpath('external/phisplit')
rmpath('external/phisplit/extern/KronPACK/src')
rmpath('external/phisplit/extern/phiks')
