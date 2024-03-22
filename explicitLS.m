function [stp,Fend,nfev,info] = explicitLS(n,m,x,d,stp,ifirst,f0,g0,Dgstp)

% Parameters

gamma = 1.999999;

% Counters

nfev = 0;

% Define ftest

ftest = Dgstp + gamma * norm(d)^2 / 2;


Findex = 1:m;
Findex(1) = ifirst;
Findex(ifirst) = 1;


% Test the sufficient descent condition (sdc) at stp = 1

sdc = true;     
for j = 1:m
    
    i = Findex(j);
		
	[f] = evalg(n,x+stp*d,i);
    nfev = nfev + 1;

	Fend(i) = f;
	
	if ( f > f0(i) + stp * ftest(i) )
		sdc = false;
		iA  = i;
		fiA = f;
		break
    end			
end

if ( sdc )
	info = 0;
	return
end

[stp,Fend,nfevbt,info] = backtrackingMO(stp,n,m,x,d,f0,g0,iA,fiA,ftest);
nfev = nfev + nfevbt;

end

%***********************************************************************
%***********************************************************************	

function [stp,Fend,nfev,info] = backtrackingMO(stp,n,m,x,d,f0,g0,iA,fiA,ftest)	

% Parameters

stpmin = 10^(-15);
sigma1 = 0.5d-1; 
sigma2 = 9.5d-1;

% Counters

outiter = 0;
nfev    = 0;
	
%-------------------------------------------------------------------
%     Main loop
%-------------------------------------------------------------------    

while (1)

	% Test the vector Armijo condition
	
	if ( outiter == 0 )
		sdc = false;
		ind = iA;
		f   = fiA;
	elseif ( infoBT == 0 )
		sdc = true;
		for i = 1:m				

			if ( i == ind ) 
                continue
            end
		
			[f] = evalg(n,x+stp*d,i);
            nfev = nfev + 1;
	        
	        Fend(i) = f;

			if ( f > f0(i) + stp * ftest(i) )
				sdc = false;
				ind = i;
				break
            end
        end
		
    end
	
	% Finish backtracking with the current point
	
	if ( sdc )
		info = 0;
		return
    end
	
	% Test if stp is too small
	
	if ( stp <= stpmin )
		stp = stpmin;
		info = 1;
		
		for i = 1:m				
			
			if ( i == ind ) 
                continue
            end
			
			[f] = evalg(n,x+stp*d,i);
            nfev = nfev + 1;
          
            Fend(i) = f;
        end
		return
    end

	outiter = outiter + 1;
	
	% Compute new trial stepsize based on f_ind
	
	while (1)
	
		% Test the descent condition for f_ind
		
		if ( f <= f0(ind) + stp * ftest(ind) )
			infoBT = 0;
			break
		end 
		
		if ( stp <= stpmin )
            infoBT = 1;
            break
        end
		
        if ( g0(ind) < 0 )
		
            stpq = ( (g0(ind) / ( (f0(ind)-f) / stp + g0(ind) ) ) / 2 ) * stp;

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
	
	Fend(ind) = f;

end

%--------------------------------------------------------------------- 
%     End of main loop
%---------------------------------------------------------------------

end