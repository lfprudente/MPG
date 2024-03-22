function [stp,Fend,nfev,info] = explicitLSone(n,x,d,stp,ind,f0,g0)

% Parameters

gamma = 1.999999;

stpmin = 10^(-15);
sigma1 = 0.5d-1; 
sigma2 = 9.5d-1;
	
% Counters

nfev = 0;

% Compute g0

%g0 = Jf * d;

% Define ftest

%[maxg0,imax] = max(g0);

% ftest = max(g0) + gamma * norm(d)^2 / 2;

ftest = g0 + gamma * norm(d)^2 / 2;  
	
    
[f] = evalg(n,x+stp*d,ind);
nfev = nfev + 1;

while (1)

    % Test the descent condition for f_ind

    if ( f <= f0 + stp * ftest )
        Fend = f;
        info = 0;
        return
    end 

    if ( stp <= stpmin )
        Fend = f;
        info = 1;
        break
    end

    if ( g0 < 0 )

        stpq = ( (g0 / ( (f0-f) / stp + g0 ) ) / 2 ) * stp;

        stpt = stpq;

        if ( stpt >= sigma1 * stp && stpt <= sigma2 * stp )
            stp = stpt;
        else
            stp = stp / 2;
        end
    else

        stp = stp / 2;

    end

    [f] = evalg(n,x+stp*d,ind);
    nfev = nfev + 1;

end


%--------------------------------------------------------------------- 
%     End of main loop
%---------------------------------------------------------------------
