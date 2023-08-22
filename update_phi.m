function phi_out = update_phi(phi,S,DT,sigma2,kappa2)
%% function for sampling phi from a truncated normal distribution 
%INPUT
% phi: current state of phi (rate of RSL changes)
% S: current state of S (rsl)
% DT: size of the temporal grids 
% sigma2: current state of sigma2
% kappa2: current state of kappa2
%OUTPUT
% phi_out: new state of phi
%%
N = numel(phi)-1;  %number of time grids
mu = zeros(N+1,1);
nu = zeros(N+1,1);
phi_out = zeros(N+1,1);
%% calculate the mean and standard deviation of phi
for i = 1:N+1
    if i == 1  %first grid point
       mu(i) = kappa2*(S(i+1)-S(i))/(DT*kappa2+sigma2);
       nu(i) = sqrt(kappa2*sigma2/(DT*kappa2+sigma2));
    elseif i > 1 && i < N+1 %internal grid points
       mu(i) = kappa2*(S(i+1)-S(i-1))/(2*DT*kappa2+sigma2);
       nu(i) = sqrt(kappa2*sigma2/(2*DT*kappa2+sigma2));
    else       %last grid point
       mu(i) = kappa2*(S(i)-S(i-1))/(DT*kappa2+sigma2);
       nu(i) = sqrt(kappa2*sigma2/(DT*kappa2+sigma2));
    end
end
%% update phi
for i = 1:N+1
    if i == 1  % first grid point
        if S(i) < S(i+1)     % rising pattern
            phi_out(i) = dtrandn_MH(phi(i),mu(i),nu(i),max(0,mu(i)-nu(i)),mu(i)+nu(i));
        elseif S(i) > S(i+1) % falling pattern
            phi_out(i) = dtrandn_MH(phi(i),mu(i),nu(i),mu(i)-nu(i),min(0,mu(i)+nu(i)));
        end 
    elseif i > 1 && i < N+1 % internal grid points
        if S(i)>S(i-1) && S(i)<S(i+1)       % monotonic rising
            phi_out(i) = dtrandn_MH(phi(i),mu(i),nu(i),max(0,mu(i)-nu(i)),mu(i)+nu(i));
        elseif S(i)<S(i-1) && S(i)>S(i+1)   % monotonic falling      
            phi_out(i) = dtrandn_MH(phi(i),mu(i),nu(i),mu(i)-nu(i),min(0,mu(i)+nu(i)));
        elseif S(i)>S(i-1) && S(i)>S(i+1)   % concave pattern
            phi_out(i) = dtrandn_MH(phi(i),mu(i),nu(i),mu(i)-nu(i),mu(i)+nu(i));  
        elseif S(i)<S(i-1) && S(i)<S(i+1)   % convex pattern
            phi_out(i) = dtrandn_MH(phi(i),mu(i),nu(i),mu(i)-nu(i),mu(i)+nu(i));  
        end
    else % last grid point
        if S(i-1) < S(i)     % rising pattern
            phi_out(i) = dtrandn_MH(phi(i),mu(i),nu(i),max(0,mu(i)-nu(i)),mu(i)+nu(i));
        elseif S(i-1) > S(i) % falling pattern
            phi_out(i) = dtrandn_MH(phi(i),mu(i),nu(i),mu(i)-nu(i),min(0,mu(i)+nu(i)));
        end    
    end
end
end