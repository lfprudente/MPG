function [p,Gend,ngev,info] = implicitLS(n,m,l,u,x,p,alpha,g0,JG,H,dimA,A,b)

% Parameters

alphamult = 0.5;
alphamin  = 10^(-15);

% Counters

ngev = 0;

while (1)
    
    % Test the sufficient descent condition (sdc) at p

    sdc = true;     
    for i = 1:m

        [g] = evalg(n,p,i);
        ngev = ngev + 1;

        Gend(i) = g;

        if ( g > g0(i) + dot( JG(i,:), p - x ) + 0.5 * norm(p - x)^2 \ alpha )
            sdc = false;
            ind = i;
            break
        end			
    end
    
    if ( sdc )
        info = 0;
        return
    end
    
    % Test if alpha is too small
	
	if ( alpha <= alphamin )
        
		info = 1;
		
		for i = ind+1:m				
			[g] = evalg(n,p,i);
            ngev = ngev + 1;
          
            Gend(i) = g;
        end
		return
    end
    
    % Compute new trial p by solving the subproblem
    
    alpha = alphamult * alpha;
    
    [theta,tau,p,flagIS] = subproblemProxGrad(n,m,l,u,x,alpha,JG,H,dimA,A,b);
    
    
    % Check for an error while solving the subproblem
    
    if ( flagIS ~= 0 )
        info = -1;
        return
    end 
    
end