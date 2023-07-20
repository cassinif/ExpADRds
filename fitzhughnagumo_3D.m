clear all
close all

addpath('integrators')
addpath('matfiles')
addpath('external/phisplit')
addpath('external/phisplit/extern/KronPACK/src')
addpath('external/phisplit/extern/phiks')

d = 3;

n = 64*ones(1,d);
a = 0*ones(1,d);
b = pi*ones(1,d);
T = 10;

deltau = 1;
deltav = 42.1887;
rho = 24.649;
a1v = 11;
a2v = 0.1;

method = 'ETD2RKds'; % ETD2RKds, ETD2RK, Lawson2b, ETD-RDP-IF, DIRK23, RK32
compute_err = true; % if true, measure error against precomputed reference solution
                     % else plot u component at time T
tol_matlab = 1e-1; % needed if using matlab ODE suite
tol_phiks = 5e-4; % needed if using ETD2RK
nsteps = 12000; % needed if not using matlab ODE suite
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

g{1} = @(t,u,v) rho*(-u.*(u.*u-1)-v);
g{2} = @(t,u,v) (rho*a1v)*(u-a2v*v);

dgdu{1}{1} = @(t,u,v) -rho*(3*(u.*u)-1); %dg1du
dgdu{1}{2} = @(t,u,v) -rho*ones(size(u)); %dg1dv
dgdu{2}{1} = @(t,u,v) (rho*a1v)*ones(size(u)); %dg2du
dgdu{2}{2} = @(t,u,v) -(rho*a1v*a2v)*ones(size(u)); %dg2dv

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
A_otimes{1}{1} = kron(speye(n(2)*n(3)),A_sp{1}{1});
A_otimes{1}{2} = kron(speye(n(3)),kron(A_sp{1}{2},speye(n(1))));
A_otimes{1}{3} = kron(A_sp{1}{3},speye(n(1)*n(2)));
A_otimes{2}{1} = kron(speye(n(2)*n(3)),A_sp{2}{1});
A_otimes{2}{2} = kron(speye(n(3)),kron(A_sp{2}{2},speye(n(1))));
A_otimes{2}{3} = kron(A_sp{2}{3},speye(n(1)*n(2)));

gvec = @(t,uvec) [g{1}(t,uvec(1:pn),uvec(pn+1:2*pn));g{2}(t,uvec(1:pn),uvec(pn+1:2*pn))];
g_if = @(u,v) gvec(NaN,[u(:);v(:)]);

load('fitzhughnagumo_3D_U0.mat')
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
    U = etd_rdp_if_3d(U0,A_otimes,g_if,nsteps,tau);
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
  load('fitzhughnagumo_3D_Uref.mat')
  normrifu = norm(Uref{1}(:),2);
  normrifv = norm(Uref{2}(:),2);
  err = norm([norm(U{1}(:)-Uref{1}(:),2)/normrifu,norm(U{2}(:)-Uref{2}(:),2)/normrifv]);
  fprintf('Error: %.3e\n',err)
else
  nn = n(1)/2;
  figure;
  surf(X{1}(:,:,nn),X{2}(:,:,nn),U{1}(:,:,nn),'edgecolor','none');
  axis equal
  view(2)
  xlabel('x_1')
  ylabel('x_2')
  colorbar
  drawnow

  idx = [2,1,3];
  lev = 0.08;
  figure;
  isosurface(permute(X{1},idx),permute(X{2},idx),X{3},permute(U{1},idx),lev)
  xlim([a(1) b(1)]); xlabel('x_1');
  ylim([a(2) b(2)]); ylabel('x_2');
  zlim([a(3) b(3)]); zlabel('x_3');
  drawnow
end

fprintf('Wall-clock time: %.2f s\n',wctime)

rmpath('integrators')
rmpath('matfiles')
rmpath('external/phisplit')
rmpath('external/phisplit/extern/KronPACK/src')
rmpath('external/phisplit/extern/phiks')
