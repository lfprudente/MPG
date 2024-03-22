function [x,info,iter,ngev,nhev,time] = ProxGrad(n,m,l,u,x,lsopt,dimA,A,b)

tic;

tol  = 10^(-4);

itermax = 200;

% Counters

ngev = 0;
nhev = 0;
iter = 0;

% Evaluate F

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

fprintf('----------------------------------------------------------------------\n')
fprintf('   Proximal Gradient Method for Composite Multiojective Optimization  \n')
fprintf('----------------------------------------------------------------------\n')
fprintf('Number of variables  : %i \n',n)
fprintf('Number of objectives : %i \n',m)
fprintf('Optimality tolerance: %.0e \n',tol)
fprintf('Line search scheme  : %i \n\n',lsopt)

% -----------
% Main Loop
% -----------

while(1)

    % Compute the Jacobian of G
    
%     if iter > 0
%         JGprev = JG;
%     end

    for i = 1:m
        JG(i,:) = evalgradg(n,x,i);
    end

    % Solve the subproblem
    
    alpha = 1;

    [theta,tau,p,flagIS] = subproblemProxGrad(n,m,l,u,x,alpha,JG,H,dimA,A,b);

    % Check for an error while solving the subproblem
    
    if ( flagIS ~= 0 )
        
        info = -1;
        
        % Stop timing 

        time = toc;

        % Print information

        fprintf('\n')
        fprintf('An error occurred while solving the subprobem H.\n')
        fprintf('Number of G_i evaluations: %i\n',ngev)
        fprintf('Number of H_i evaluations: %i\n',nhev)
        fprintf('CPU time(s)                   : %.1f \n',time)
        return
    end 

    % Print information

    if ( mod(iter,10) == 0 )
        fprintf('\n')
        fprintf('%-5s   %-8s       %-8s   %-8s %-8s %-8s \n','it','|theta|','fun1','fun2','IS','LS')
    end
    if ( iter == 0 )
        fprintf('%5d   %8.2e    %+8.2e   %+8.2e   %-8i %-8s\n',iter,abs(theta),F(1),F(2),flagIS,'-')
    else
        fprintf('%5d   %8.2e    %+8.2e   %+8.2e   %-8i %-8i\n',iter,abs(theta),F(1),F(2),flagIS,flagLS)
    end
    
    % -----------
    % Stopping criteria
    % -----------

    % Test optimality 

    if ( abs(theta) <= tol )
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
    % Iterate
    % -----------

    iter = iter + 1;

    % Define the search direction

    if ( lsopt == 1 || lsopt == 2 ), d = p - x;, end
    
    % Explicit Multiobjective Proximal Gradient

    if ( lsopt == 1 )
        
        stp = 1;

        DG0 = JG * d;

        [maxDG0,imax] = max(DG0);

        [stp,Gtrialimax,ngevLS,flagLS] = explicitLSone(n,x,d,stp,imax,G(imax),DG0(imax));
        ngev = ngev + ngevLS;
        
        xtrial = x + stp * d;
        
        Ftrial(1:m) = 0;
        Gtrial(1:m) = 0;
        Htrial(1:m) = 0;
        
        Gtrial(imax) = Gtrialimax;
        
        Fdec = true;
        for i = 1:m
            
            if ( i == imax ) 
                continue
            end
            
            [Gtrial(i)] = evalg(n,xtrial,i);
            ngev = ngev + 1;
            
            [Htrial(i),flag] = evalh(n,xtrial,i,A,b);
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
            
            Ftrial(i) = Gtrial(i) + Htrial(i);
            if ( Ftrial(i) > F(i) )
                Fdec = false;
                ifirst = i;
                break;
            end
            
        end
        
        if ( Fdec )
            
            [Htrial(imax),flag] = evalh(n,xtrial,imax,A,b);
            nhev = nhev + 1;
            
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
            
            Ftrial(imax) = Gtrial(imax) + Htrial(imax);
            
            % Update x
            
            x = xtrial;
            G = Gtrial;
            H = Htrial;
            F = Ftrial;
        else

            [stp,G,ngevLS,flagLS] = explicitLS(n,m,x,d,stp,ifirst,G,DG0,DG0);
            ngev = ngev + ngevLS;
            
            % Check for an error in the line seach procedure
    
            if ( flagLS == -1 )

                info = -1;

                % Stop timing 

                time = toc;

                % Print information

                fprintf('\n')
                fprintf('An error occurred in the line seach procedure.\n')
                fprintf('Number of G_i evaluations: %i\n',ngev)
                fprintf('Number of H_i evaluations: %i\n',nhev)
                fprintf('CPU time(s)                   : %.1f \n',time)
                return
            end 
            
            % Update x
            
            x = x + stp * d;
            
            for i = 1:m
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
            
            % Evaluate F

            F = G + H;

        end

    end
    
    % Multiobjective Proximal Gradient with Armijo line search
    
    if ( lsopt == 2 )

        [stp,G,H,ngevLS,nhevLS,flagLS] = armijo(n,m,x,d,F,JG,tau,A,b);
        ngev = ngev + ngevLS;
        nhev = nhev + nhevLS;
        
         % Check for an error in the line seach procedure

        if ( flagLS == -1 )

            info = -1;

            % Stop timing 

            time = toc;

            % Print information

            fprintf('\n')
            fprintf('An error occurred in the line seach procedure.\n')
            fprintf('Number of G_i evaluations: %i\n',ngev)
            fprintf('Number of H_i evaluations: %i\n',nhev)
            fprintf('CPU time(s)                   : %.1f \n',time)
            return
        end 
        
        % Update x
        
        x = x + stp * d;
        
        % Evaluate F
        
        F = G + H;
       
    end
    
    % Implicit Multiobjective Proximal Gradient
    
    if ( lsopt == 3 )

        [p,G,ngevLS,flagLS] = implicitLS(n,m,l,u,x,p,alpha,G,JG,H,dimA,A,b);
        ngev = ngev + ngevLS;
        
         % Check for an error in the line seach procedure
    
        if ( flagLS == -1 )

            info = -1;

            % Stop timing 

            time = toc;

            % Print information

            fprintf('\n')
            fprintf('An error occurred in the line seach procedure.\n')
            fprintf('Number of G_i evaluations: %i\n',ngev)
            fprintf('Number of H_i evaluations: %i\n',nhev)
            fprintf('CPU time(s)                   : %.1f \n',time)
            return
        end 
        
        % Update x
        
        x = p;
        
        % Evaluate H
         
        for i = 1:m
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
            
            % Evaluate F
            
            F = G + H;
            
        end

    end
    
end

% -----------
% End of Main Loop
% -----------
