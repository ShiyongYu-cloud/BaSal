function R_hat = convergence(mcmcsamples,chain_name)
%% function for monitoring convergence of Markov Chains
% This is an implementation of Andrew Gelman's algorithm in chapter for Gilks, Richardson, and
% Spiegelhalter book
%INPUT   
% mcmcsample: nsamples*K*nchains array containing results from MCMC run
% chain_name: name of the chain
%OUTPUT
% R_hat: potential scale reduction
%%
[nsamples,K,nchains] = size(mcmcsamples);
r = zeros(1,K);
X = zeros(nsamples,nchains);
for i = 1:K             % loop over the parameters 
    for j = 1:nchains   % loop over the chains
        X(:,j) = mcmcsamples(:,i,j);  % extract simulations for each parameter
    end    
    B = var(mean(X),0);  % calculate the between-sequence variance
    W = mean(var(X,0));  % calculate the within-sequence variance
    V = (nsamples-1)*W/nsamples + B/nsamples; % conservative estimate of variance
    r(i) = V/W; % calcualte potential scale reduction for each parameter
end
r(isnan(r)) = []; % remove NaN
R_hat = power(prod(r),1/K); % calculate the geometric mean
if R_hat <= 1.2
    disp(['R_hat = ' num2str(R_hat) ' and the chain of ' chain_name ' converged!']);
else
    disp(['R_hat = ' num2str(R_hat) ' and the chain of ' chain_name ' not converged!']);
end
end