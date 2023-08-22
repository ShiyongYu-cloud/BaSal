function kappa2_out = update_kappa2(phi)
%% function for sampling sigma2 from the inverse Gamma distribution
%INPUT:
% phi: current state of phi
%OUTPUT
% kappa2_out: new state of kappa2
%%
N = numel(phi)-1; %number of time grids
alpha = (N+1)/2;
beta = sum(phi.*phi)/2; 
kappa2_out = 1/gamrnd(alpha,1/beta,1,1);
end