function [x,info,iter,ngev,nhev,time] = ProxGradAcc(n,m,l,u,x,lsopt,dimA,A,b)

tic;

% Initial parameters

tol  = 10^(-4);

alpha     = 1;
alphamult = 0.5;

itermax      = 200;
backtrackmax = 100;

x = min( u, max(x,l) );

% Set the initial momentum term

y = x;
t = 1;

% Counters

ngev = 0;
nhev = 0;
iter = 0;
backtrack = 0;

% Evaluate F at x

for i = 1:m
    G(i) = evalg(n,x,i);
    ngev = ngev + 1;

    [H(i),flag] = evalh(n,x,i,A,b);
    nhev = nhev + 1;
    
    % Check for an error while evaluating H
    
    if ( flag ~= 0 )
        
        info = -1;
        
        % Stop timing 

        time = toc;

        % Print information

        fprintf('\n')
        fprintf('An error occurred while evaluating H.\n')
        fprintf('Number of G_i evaluations: %i\n',ngev)
        fprintf('Number of H_i evaluations: %i\n',nhev)
        fprintf('CPU time(s)                   : %.1f \n',time)
        return
    end 
end

F = G + H;

fprintf('----------------------------------------------------------------------------------\n')
fprintf('   Accelerated Proximal Gradient Method for Composite Multiojective Optimization  \n')
fprintf('----------------------------------------------------------------------------------\n')
if ( lsopt == 1 )
    fprintf('Version              : Accelerated \n')
elseif ( lsopt == 2 )
    fprintf('Version              : Normal \n')
end
fprintf('Number of variables  : %i \n',n)
fprintf('Number of objectives : %i \n',m)
fprintf('Optimality tolerance: %.0e \n',tol)

fprintf('\n')
fprintf('%-5s     %-5s     %-8s       %-8s   %-8s %-8s  %-8s\n','it','|x-y|','|theta|','fun1','fun2','IS','L')
fprintf('%5d       %-4s       %-8s %+8.2e   %+8.2e   %-8s %-8s\n',iter,'-','-',F(1),F(2),'-','-')

% -----------
% Main Loop
% -----------

while(1)
    
    % -----------
    % Iterate
    % -----------

    iter = iter + 1;

    if ( lsopt == 1 )
    
        % Evaluate G at y 
        
        for i = 1:m
            Gy(i) = evalg(n,y,i);
            ngev = ngev + 1;
        end

        FxGy = F - Gy;
    elseif ( lsopt == 2 )

        FxGy = H;
    end
    
    % Compute the Jacobian of G at y

    for i = 1:m
        JG(i,:) = evalgradg(n,y,i);
    end

    % Salve the current value of x and F

    xprev = x;
    Fprev = F;

    % Backtracking scheme to estimate the Lipschitz constant

    while(1)

        % Solve the subproblem

        [theta,~,p,flagIS] = subproblemProxGrad(n,m,l,u,y,alpha,JG,FxGy,dimA,A,b);
    
        % Check for an error while solving the subproblem
        
        if ( flagIS ~= 0 )
            
            info = -2;
            
            % Stop timing 
    
            time = toc;
    
            % Print information
    
            fprintf('\n')
            fprintf('An error occurred while solving the subproblem.\n')
            fprintf('Number of G_i evaluations: %i\n',ngev)
            fprintf('Number of H_i evaluations: %i\n',nhev)
            fprintf('CPU time(s)                   : %.1f \n',time)
            return
        end 

        x = p;

        % Evaluate F at x
        
        for i = 1:m
            G(i) = evalg(n,x,i);
            ngev = ngev + 1;
        
            [H(i),flag] = evalh(n,x,i,A,b);
            nhev = nhev + 1;
            
            % Check for an error while evaluating H
            
            if ( flag ~= 0 )
                
                info = -1;
                
                % Stop timing 
        
                time = toc;
        
                % Print information
        
                fprintf('\n')
                fprintf('An error occurred while evaluating H.\n')
                fprintf('Number of G_i evaluations: %i\n',ngev)
                fprintf('Number of H_i evaluations: %i\n',nhev)
                fprintf('CPU time(s)                   : %.1f \n',time)
                return
            end 
        end
        
        F = G + H;

        if ( max(F - Fprev) <= theta + 10^(-12) || norm(x-y,inf)/max(1,norm(y,inf)) <= tol )
            break
        else
            backtrack = backtrack + 1;
            alpha = alphamult * alpha;

            if ( backtrack > backtrackmax )
                
                info = -3;
                
                % Stop timing 
        
                time = toc;
        
                % Print information
        
                fprintf('\n')
                fprintf('Backtracking failed to find a suitable stepsize.\n')
                fprintf('Number of G_i evaluations: %i\n',ngev)
                fprintf('Number of H_i evaluations: %i\n',nhev)
                fprintf('CPU time(s)                   : %.1f \n',time)
                return
            end
        end
    end
    
    % Compute norm(x-y,inf)/max(1,norm(y,inf)

    supnormxy = norm(x-y,inf)/max(1,norm(y,inf));

    % Print information

    if ( iter > 0 && mod(iter,10) == 0 )
        fprintf('\n')
        fprintf('%-5s     %-5s     %-8s       %-8s   %-8s %-8s  %-8s\n','it','|x-y|','|theta|','fun1','fun2','IS','L')
    end
    fprintf('%5d    %8.2e   %8.2e    %+8.2e   %+8.2e   %-8i %8.2e\n',iter,supnormxy,abs(theta),F(1),F(2),flagIS,1/alpha)

    % -----------
    % Stopping criteria
    % -----------

    % Test optimality 

     if ( supnormxy <= tol )
        info = 0;

        % Stop timing 

        time = toc;

        % Print information

        fprintf('\n')
        fprintf('Solution was found.\n')
        fprintf('Number of G_i evaluations: %i\n',ngev)
        fprintf('Number of H_i evaluations: %i\n',nhev)
        fprintf('CPU time(s)                   : %.1f \n',time)

        return
     end
            
    % Test the maximum number of iterations

    if ( iter >= itermax )
        info = 1;

        % Stop timing 

        time = toc;

        % Print information

        fprintf('\n')
        fprintf('The number of maximum iterations was reached.\n')
        fprintf('CPU time(s): %.1f \n',time)

        return
    end 

    % -----------
    % Prepare the next iterate
    % -----------

    if ( lsopt == 1 )

        % Update the initial momentum term

        tprev = t;
        t = sqrt( tprev^2 + 0.25 ) + 0.5;
        gammak = ( tprev - 1 ) / t;
        y = x + gammak * ( x - xprev );

    elseif ( lsopt == 2 )

        y = x;
    end

end

% -----------
% End of Main Loop
% -----------