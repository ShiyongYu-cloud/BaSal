function [S,phi,sigma2,kappa2] = Gibbs(S,phi,T,G,DT,Y,age_scale)
% function for updating other model parameters using the Gibbs samplers
%INPUT:
% S:      current state of S (rsl)
% phi:    current state of phi (rate of rsl changes)
% T:      modeled calendar ages of the SLIPs
% G;      temporal grid points
% DT:     size of the time grid  
% Y:      structure containing altitudinal info of the SLIPs
% age_sale: scale of year (BC/AD,BP,B2K)
%OUTPUT:
% S:  new state of S   
% phi: new state of phi
% sigma2: new state of sigma2
% kappa2: new state of kappa2
%% 
sigma2 = update_sigma2(S,phi,DT);
kappa2 = update_kappa2(phi);        
phi = update_phi(phi,S,DT,sigma2,kappa2);
S = update_S(S,phi,sigma2,T,G,DT,Y,age_scale);
end