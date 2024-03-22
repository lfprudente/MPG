function [dimA,A,b] = datas(n,m,delta)

dimA(1:m) = 2*n;

l = -10;
u =  10;

for ind = 1:m
    B = l + rand(n) * (u-l);
    A{ind} = [B; -B];
    b{ind}   = delta * ones(2*n,1);
end