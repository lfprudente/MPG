function [h,info] = evalh(n,x,ind,A,b)

options.Algorithm = 'dual-simplex';
options.Display = 'off';
options.OptimalityTolerance = 1e-10;

[xopt,fmin,flag] = linprog(-x,A{ind},b{ind},[],[],-inf(n,1),inf(n,1),options);

h = -fmin;

info = -1;
if ( flag == 1 )
    info = 0;
    return
end