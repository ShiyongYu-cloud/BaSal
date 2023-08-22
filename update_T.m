function T_out = update_T(T,X,A,B,age_scale,CalCurves)
%% function for updating the calendar age using the Metropolis-Hastings algorithm
%INPUT
% T: initial state of the modeled calendar ages
% X: structure containing information about ages of the SLIPs
% A: early boundary of the age sequence
% B: late boundary of the age sequence
% age_scale: scale of year to be reported (BC/AD, BP, B2K)
% CalCurve: sturcture containing the calibration curves for 14C ages in BP
%OUTPUT
% T_out: updated state of the modeled calendar ages
%%
% propose a random-walk move of the calendar ages
N = length(T);
T_star = round(T+5*(2*rand(N,1)-1));
true = length(unique(T_star)) == length(T_star);% ensure unique ages 
while(~true)
    T_star = round(T+5*(2*rand(N,1)-1));  
    true = length(unique(T_star)) == length(T_star);
end
% calculate the loglikelihood
llk_num = age_like(X,T_star,A,B,age_scale,CalCurves);
llk_den = age_like(X,T,A,B,age_scale,CalCurves);
% calculate the log prior
prior_num = age_prior(T_star,A,B,age_scale);
prior_den = age_prior(T,A,B,age_scale);
% calculate the log posterior ratio
omega = (llk_num+prior_num)-(llk_den+prior_den);
omega = min(omega,0);
rho = log(rand(1,1));
if omega > rho
   T_out = T_star;
else
   T_out = T;
end   
end