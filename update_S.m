function S_out = update_S(S,phi,sigma2,T,G,DT,Y,age_scale)
%% function for sampling S from a truncated normal distribution 
%INPUT:
% S:      current state of S (rsl)
% phi:    current state of phi (rate of rsl changes)
% sigma2: currant state of sigma2 (stepsize of the Brownian motion)
% T:      modeled calendar ages of the SLIPs
% G;      temporal grid points
% DT:     size of the time grid  
% Y:      structure containing altitudinal info of the SLIPs
% age_sale: scale of year (BC/AD,BP,B2K)
%OUTPUT:
% S_out:  new state of S   
%%
N = numel(G)-1; % number of time grids
y = zeros(N+1,1);
v2 = zeros(N+1,1);
mv = zeros(N+1,1);
nv = zeros(N+1,1);
S_out = zeros(N+1,1);
%% round T to the nearest time grid points in G
[T_G,id,~] = unique(round(T/DT)*DT,'stable');
y_G = Y.rsl(id);
v2_G = Y.rsl_err(id).^2;
%% read observational data to the nearest time grid points 
for i = 1:N+1   %loop over all time grid points
    for j = 1:numel(T_G)
        if G(i) == T_G(j) 
           y(i) = y_G(j);
           v2(i) = v2_G(j);
           break;
        end
    end
end 
%% calculate the mean and standard deviation of S
for i = 1:N+1
    if i == 1                      %first grid point
       mv(i) = S(i+1) - phi(i+1)*DT;
       nv(i) = sqrt(DT*sigma2);
    elseif i > 1 && i < N+1        %internal grid points
       if ismember(G(i),T_G) == 1    %grid points with rsl data
          mv(i) = (DT*sigma2*y(i)+v2(i)*(S(i-1)+phi(i-1)*DT+S(i+1)-phi(i+1)*DT))/(DT*sigma2+2*v2(i));
          nv(i) = sqrt(DT*sigma2*v2(i)/(DT*sigma2+2*v2(i)));
       else                        %internal grid points without rsl data
          mv(i) = (S(i-1)+phi(i-1)*DT+S(i+1)-phi(i+1)*DT)/2;
          nv(i) = sqrt(DT*sigma2/2);
       end
    else                           %last grid point 
       mv(i) = S(i-1) + phi(i-1)*DT;
       nv(i) = sqrt(DT*sigma2);
     end
end
%% update S
for i = 1:N+1
    if i == 1  % first grid point
        if S(i) < S(i+1)     % rising pattern
            S_out(i) = dtrandn_MH(S(i),mv(i),nv(i),S(i)-phi(i)*DT,S(i+1));
        elseif S(i) > S(i+1) % falling pattern
            S_out(i) = dtrandn_MH(S(i),mv(i),nv(i),S(i+1),S(i)-phi(i)*DT);
        end 
    elseif i > 1 && i < N+1  % internal grid points
        if S(i)>S(i-1) && S(i)<S(i+1)       % monotonic rising
            S_out(i) = dtrandn_MH(S(i),mv(i),nv(i),S_out(i-1),S(i+1));
        elseif S(i)<S(i-1) && S(i)>S(i+1)   % monotonic falling      
            S_out(i) = dtrandn_MH(S(i),mv(i),nv(i),S(i+1),S_out(i-1));
        elseif S(i)>S(i-1) && S(i)>S(i+1)   % concave pattern
            S_out(i) = dtrandn_MH(S(i),mv(i),nv(i),max(S_out(i-1),S(i+1)),min(S(i)+phi(i-1)*DT,S(i)-phi(i+1)*DT));  
        elseif S(i)<S(i-1) && S(i)<S(i+1)   % convex pattern
            S_out(i) = dtrandn_MH(S(i),mv(i),nv(i),max(S(i)+phi(i-1)*DT,S(i)-phi(i+1)*DT),min(S_out(i-1),S(i+1)));  
        end
    else % last grid point
        if S(i-1) < S(i)     % rising pattern
            S_out(i) = dtrandn_MH(S(i),mv(i),nv(i),S_out(i-1),min(S(i)+phi(i)*DT,0));
        elseif S(i-1) > S(i) % falling pattern
            S_out(i) = dtrandn_MH(S(i),mv(i),nv(i),max(S(i)+phi(i)*DT,0),S_out(i-1));
        end    
    end
end
% ensure that RSL at the origin is zero if B is the origin
if strcmpi(age_scale,'BC/AD') == 1
   if G(end) == 2000
       S_out(end) = 0;
   end
else
   if G(end) == 0
       S_out(end) = 0;  
   end  
end
end