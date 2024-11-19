function [theta,tau,p,flag] = subproblemProxGrad(n,m,l,u,x,alpha,JG,H,dimA,A,b)

flag = -1;

sdimA = sum(dimA);

Hess = zeros(1+n+sdimA);
Hess(2:n+1,2:n+1) = eye(n) / alpha;

c = zeros(1+n+sdimA,1);
c(1) = 1;
c(2:n+1) = - x / alpha;

index(1) = n + 1;
for ind = 2:m+1
    index(ind) = index(ind-1) + dimA(ind-1);
end

Asubprob   = zeros(m,1+n+sdimA);
Aeqsubprob = zeros(n*m,1+n+sdimA);

for ind = 1:m
    Asubprob(ind,1:1+n) = [-1 JG(ind,:)];
    Asubprob(ind,index(ind)+1:index(ind+1)) = b{ind};
    Aeqsubprob((ind-1)*n+1:ind*n,2:n+1) = -eye(n);
    Aeqsubprob((ind-1)*n+1:ind*n,index(ind)+1:index(ind+1)) = A{ind}';
end

for ind = 1:m
    bsubprob(ind) = H(ind) + dot(JG(ind,:),x); 
end
bsubprob = bsubprob';

beqsubprob(1:n*m) = 0;
beqsubprob = beqsubprob';

lb = [-inf; l; zeros(sdimA,1)];
ub = [ inf; u; inf(sdimA,1)  ];

%options.Algorithm = 'interior-point-convex';
options.Display = 'off';
options.OptimalityTolerance = 1e-10;

try
    [xopt, fmin, info] = quadprog(Hess, c, Asubprob, bsubprob, Aeqsubprob, beqsubprob,lb, ub, [0; x],options);

    if ( info == 1 )
        flag = 0;
    end
    
    p = xopt(2:n+1);
    tau   = xopt(1);
    theta = xopt(1) + 0.5 * norm( p - x )^2 / alpha;

catch exception

    flag = -1;

    theta = NaN;
    tau = NaN;
    p = NaN;

end

