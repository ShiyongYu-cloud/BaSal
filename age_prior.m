function p = age_prior(T,A,B,age_scale)
%% function for calculating the prior probability of the modeled calendar ages
%INPUT
% T: modeled calendar ages of the SLIPs
% A: early boundary of the age sequence
% B: late boundary of the age sequence
% age_scale: scale of year (BC/AD, BP, B2K)
%OUTPUT
% p: prior probability of the modeled calendar ages assuming uniform distribution on [A, B] 
%%
N = length(T);
if strcmpi(age_scale,'BC/AD') == 1  
    if all(A < T) && all(T < B)
        f = power(B-A,-N);
    else 
        f = 0;
    end    
    if all(diff(T)>0)
        I = 1;
    else
        I = 0;
    end
else 
    if all(B < T) && all(T < A)
        f = power(A-B,-N);
    else 
        f = 0;
    end    
    if all(diff(T)<0)
        I = 1;
    else
        I = 0;
    end
end
p = log(f) + log(I);
end