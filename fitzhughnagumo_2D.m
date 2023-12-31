clear all
close all

addpath('integrators')
addpath('matfiles')
addpath('external/phisplit')
addpath('external/phisplit/extern/KronPACK/src')
addpath('external/phisplit/extern/phiks')

d = 2;

n = 100*ones(1,d);
a = 0*ones(1,d);
b = pi*ones(1,d);
T = 10;

deltau = 1;
deltav = 42.1887;
rho = 65.731;
a1v = 11;
a2v = 0.1;

method = 'ETD2RKds'; % ETD2RKds, ETD2RK, Lawson2b, ETD-RDP-IF, DIRK23, RK32
compute_err = true; % if true, measure error against precomputed reference solution
                    % else plot u component at time T
compute_indic = false; % if true compute indicators (relevant if compute_err=false
                       % and method='ETD2RKds')
tol_matlab = 5e-7; % needed if using matlab ODE suite
tol_phiks = 1e-3; % needed if using ETD2RK
nsteps = 20000; % needed if not using matlab ODE suite
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
A_otimes{1}{1} = kron(speye(n(2)),A_sp{1}{1});
A_otimes{1}{2} = kron(A_sp{1}{2},speye(n(1)));
A_otimes{2}{1} = kron(speye(n(2)),A_sp{2}{1});
A_otimes{2}{2} = kron(A_sp{2}{2},speye(n(1)));

gvec = @(t,uvec) [g{1}(t,uvec(1:pn),uvec(pn+1:2*pn));g{2}(t,uvec(1:pn),uvec(pn+1:2*pn))];
g_if = @(u,v) gvec(NaN,[u(:);v(:)]);

load('fitzhughnagumo_2D_U0.mat')
u0 = [U0{1}(:);U0{2}(:)];

fprintf('Method: %s\n',method)
switch method
  case 'ETD2RKds'
    if compute_indic
      tic
      [U,Umean,Uinc] = etd2rkds(U0,A,F,g,nsteps,tau,x);
      wctime = toc;
    else
      tic
      U = etd2rkds(U0,A,F,g,nsteps,tau);
      wctime = toc;
    end
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
  load('fitzhughnagumo_2D_Uref.mat')
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

if compute_indic
  trange = 0:tau:T;
  figure;
  plot(trange,Umean,'r-')
  xlabel('t')
  ylabel('<U_n>')
  title('Spatial mean')
  legend(sprintf('tau=%.2e',tau))
  figure;
  semilogy(trange(2:end),Uinc,'r-')
  xlabel('t')
  ylabel('||U_{n+1}-U_n||_F')
  title('Time increment')
  legend(sprintf('tau=%.2e',tau))
end

fprintf('Wall-clock time: %.2f s\n',wctime)

rmpath('integrators')
rmpath('matfiles')
rmpath('external/phisplit')
rmpath('external/phisplit/extern/KronPACK/src')
rmpath('external/phisplit/extern/phiks')
