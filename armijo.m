function [stp,Gend,Hend,ngev,nhev,info] = armijo(n,m,x,d,f0,Jf,theta,A,b)

% Parameters

ftol = 10^(-4);

% Counters

ngev = 0;
nhev = 0;

% Compute g0

g0 = Jf * d;

% Define ftest

ftest = - ftol * abs(theta);	

% Test the sufficient descent condition (sdc) at stp = 1

stp = 1;

sdc = true;     
for i = 1:m
    
    % Evaluate g
		
	[g] = evalg(n,x+stp*d,i);
    ngev = ngev + 1;
    
    % Evaluate h

    [h,flag] = evalh(n,x+stp*d,i,A,b);
    nhev = nhev + 1;
    
     % Check for an error while evaluating h

    if ( flag ~= 0 )
        info = -1;
        return
    end 
    
    f = g + h;
	
	Gend(i) = g;
    Hend(i) = h;
	
	if ( f > f0(i) + ftest * stp )
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

[stp,Gend,Hend,ngevbt,nhevbt,info] = backtrackingMO(stp,n,m,x,d,f0,g0,iA,fiA,ftest,theta,A,b);
ngev = ngev + ngevbt;
nhev = nhev + nhevbt;

end

%***********************************************************************
%***********************************************************************	

function [stp,Gend,Hend,ngev,nhev,info] = backtrackingMO(stp,n,m,x,d,f0,g0,iA,fiA,ftest,theta,A,b)	

% Parameters

stpmin = 10^(-15);
sigma1 = 0.5d-1; 
sigma2 = 9.5d-1;

% Counters

outiter = 0;
ngev    = 0;
nhev    = 0;
	
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
		
			[g] = evalg(n,x+stp*d,i);
            ngev = ngev + 1;

            [h,flag] = evalh(n,x+stp*d,i,A,b);
            nhev = nhev + 1;
            
            % Check for an error while evaluating h

            if ( flag ~= 0 )
                info = -1;
                return
            end 
            
            f = g + h;
	        
	        Gend(i) = g;
            Hend(i) = h;

			if ( f > f0(i) + ftest * stp )
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
			
			[g] = evalg(n,x+stp*d,i);
            ngev = ngev + 1;

            [h,flag] = evalh(n,x+stp*d,i,A,b);
            nhev = nhev + 1;
            
            % Check for an error while evaluating h

            if ( flag ~= 0 )
                info = -1;
                return
            end 
            
            f = g + h;
            
            Gend(i) = g;
            Hend(i) = h;
        end
		return
    end

	outiter = outiter + 1;
	
	% Compute new trial stepsize based on f_ind
	
	while (1)
	
		% Test Armijo condition for f_ind
		
		if ( f <= f0(ind) + ftest * stp )
			infoBT = 0;
			break
		end 
		
		if ( stp <= stpmin )
            infoBT = 1;
            break
        end
		
% 		stpq = ( (g0(ind) / ( (f0(ind)-f) / stp + g0(ind) ) ) / 2 ) * stp;
% 		
% 		stpt = stpq;
% 		
% 		if ( stpt >= sigma1 * stp && stpt <= sigma2 * stp )
% 			stp = stpt;
% 		else
% 			stp = stp / 2;
%         end

        stp = stp / 2;
		
		[g] = evalg(n,x+stp*d,ind);
        ngev = ngev + 1;

        [h,flag] = evalh(n,x+stp*d,ind,A,b);
        nhev = nhev + 1;
        
        % Check for an error while evaluating h

        if ( flag ~= 0 )
            info = -1;
            return
        end 
        
        f = g + h;
    end
	
	Gend(ind) = g;
    Hend(ind) = h;

end

%--------------------------------------------------------------------- 
%     End of main loop
%---------------------------------------------------------------------

end