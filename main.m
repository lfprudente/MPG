clear
clc

global problem

rng(123);

% Choose the problem to be solved

problem = 'AP1';

% Initial parameters

[n,m,l,u,x0] = inip;                                                                                                                                                            

% Define the uncertainty parameter

delta = (8*rand+2)/100 * norm(x0);

% Set problem data

[dimA,A,b] = datas(n,m,delta);

% Choose de line search scheme to be used

% AlgOpt = 1: Explicit Multiobjective Proximal Gradient
% AlgOpt = 2: Multiobjective Proximal Gradient with Armijo line search
% AlgOpt = 3: Implicit Multiobjective Proximal Gradient
% AlgOpt = 4: Accelerated Multiobjective Proximal Gradient
% AlgOpt = 5: Normal Multiobjective Proximal Gradient

AlgOpt = 1;

% Call the solvers

if ( AlgOpt == 1 ) [x,info] = ProxGrad(n,m,l,u,x0,1,dimA,A,b);, end

if ( AlgOpt == 2 ) [x,info] = ProxGrad(n,m,l,u,x0,2,dimA,A,b);,end

if ( AlgOpt == 3 ) [x,info] = ProxGrad(n,m,l,u,x0,3,dimA,A,b);,end

if ( AlgOpt == 4 ), [x,info,iter,ngev,nhev,time] = ProxGradAcc(n,m,l,u,x0,1,dimA,A,b);, end

if ( AlgOpt == 5 ), [x,info,iter,ngev,nhev,time] = ProxGradAcc(n,m,l,u,x0,2,dimA,A,b);, end