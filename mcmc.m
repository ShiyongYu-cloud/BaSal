function [T_samples,S_samples,phi_samples,sigma2_samples,kappa2_samples] = mcmc(T,X,Y,age_scale,CalCurves,G,nsamples,burn_in,thin)
%% function for sampling the posterior distribution of T, S, and phi using MCMC method
%INPUT
% T: initial values of the calendar ages of the SLIPs
% S: initial values of RSLs 
% phi: initial values of rate of RSL changes
% sigma2: stepsize of the Brownnian motion
% kappa2: vriance of phi prior
% X: structure containing the age info of SLIPs
% Y: structure containing the altitudinal info of SLIPs
% age_scale: scale of year (BC/AD, BP, B2K)
% CalCurves: structure containing the extraced calibration curves
% G: temporal grid points
% nsamples: number of sample to be kept during MCMC run
% burn_in: number of steps from which samples will be kept during MCMC run
% thin: steps every of which samples will be kept during MCMC run
%OUTPUT
% T_samples: calendar ages generated from MCMC run
% S_samplesL RSLs generated from MCMC run
% phi_samples: rate of RSL changes generated from MCMC run
% sigma2_samples: MCMC samples of sigma2
% kappa2_samples: MCMC samples of kappa2
%%
T_samples = zeros(nsamples,length(T));
S_samples = zeros(nsamples,length(G));
phi_samples = zeros(nsamples,length(G));
sigma2_samples = zeros(nsamples,1);
kappa2_samples = zeros(nsamples,1);
A = G(1);
B = G(end);
DT = unique(abs(diff(G)));
for i = 1-burn_in:nsamples*thin
    % update the calendar ages of the SLIPs
    T = update_T(T,X,A,B,age_scale,CalCurves);
    if (i > 0) && (mod(i,thin) == 0)
        % update other model parameters
        [S,phi] = initials(T,Y,age_scale,G);
        [S,phi,sigma2,kappa2] = Gibbs(S,phi,T,G,DT,Y,age_scale);
        % save the results
        T_samples(i/thin,:) = T'; 
        S_samples(i/thin,:) = S';
        phi_samples(i/thin,:) = phi';
        sigma2_samples(i/thin) = sigma2;
        kappa2_samples(i/thin) = kappa2;
    end   
end
end