clear
clc

global problem

rng(123);

% Choose the problem to be solved

problem = 'BK1';

% Initial parameters

[n,m,l,u,x0] = inip;                                                                                                                                                            

% Define the uncertainty parameter

delta = (8*rand+2)/100 * norm(x0);

% Set problem data

[dimA,A,b] = datas(n,m,delta);

% Choose de line search scheme to be used

% lsopt = 1: Explicit Multiobjective Proximal Gradient
% lsopt = 2: Multiobjective Proximal Gradient with Armijo line search
% lsopt = 3: Implicit Multiobjective Proximal Gradient

lsopt = 1;

% Call the solvers

if ( lsopt == 1 ) [x,info] = ProxGrad(n,m,l,u,x0,lsopt,dimA,A,b);, end

if ( lsopt == 2 ) [x,info] = ProxGrad(n,m,l,u,x0,lsopt,dimA,A,b);,end

if ( lsopt == 3 ) [x,info] = ProxGrad(n,m,l,u,x0,lsopt,dimA,A,b);,end