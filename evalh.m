function [h,info] = evalh(n,x,ind,A,b)

options.Algorithm = 'dual-simplex-legacy';
options.Display = 'off';
options.OptimalityTolerance = 1e-12;

[xopt,fmin,flag] = linprog(-x,A{ind},b{ind},[],[],-inf(n,1),inf(n,1),options);

if ( flag == 1 )
    info = 0;
    h = -fmin;
else
    flag
    
    info = -1;
    h = NaN;;
end