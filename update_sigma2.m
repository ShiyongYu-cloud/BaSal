function sigma2_out = update_sigma2(S,phi,DT)
%% function for sampling sigma2 from the inverse Gamma distribution
%INPUT:
% S: current state of S
% phi: current state of phi
% DT: size of the model time grid 
%OUTPUT
% sigma2_out: new state of sigma2
%%
N = numel(S)-1; %number of time grids
alpha = (N+1)/2;
beta = 0;
for i = 1:N+1
    if i == 1
        beta = beta + (S(i)-(S(i+1)-phi(i+1)*DT))^2/DT;
    elseif i > 1 && i < N+1
        beta = beta + (S(i)-(S(i-1)+phi(i-1)*DT+S(i+1)-phi(i+1)*DT)/2)^2/DT;
    else
        beta = beta + (S(i)-(S(i-1)+phi(i-1)*DT))^2/DT;
    end
end    
sigma2_out = 1/gamrnd(alpha,1/beta,1,1);
end
