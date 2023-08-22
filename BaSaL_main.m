%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BaSaL: A MATLAB package for inferring the pattern and rate of past     % 
%              RSL changes                                               %
% Copyright: Shiyong Yu                                                  %  
% E-mail: shiyong.yu@gmail.com                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clc;
clear;
%% set up model parameters
data_file = 'SLIPS_template.txt';  % specify data file name
age_scale = 'BP';                 % specify age scale (BC/AD,BP,B2K)
DT = 50;                           % specify the size of time grid
delta = 10;        % specify the nearest years the modeled calendar ages to be rounded  
alpha = 0.05;     % specify the significance level for estimating confidence interval
nsamples = 1000;   % specify the number of samples to be kept for each chain
nchains = 3;      % specify the number of chains to be run for convergence diagnosis
burn_in = 1000;   % specify the burn-in period of the Markov chains
thin = 20;        % specify the thinning interval of the Markov chains
%% read and preprocess data (convert to F14C for 14C ages and to the specified age scale for non-14C ages)
[X,Y] = read_data(data_file);             %read data  
X = data_process(X,age_scale);            %preprocess the data   
[A,B] = age_bound(X,age_scale);           %estimate age boundaries of the model time domain 
B = 0; %reset B = 0 on the BP scale for data containling 14C ages
%B = 1950; %reset B = 1950 on the BC/AD scale
%B = 0;    %reset B = 0 on the B2K scale for data that do not contain 14C ages; otherwise, B = 50
G = make_grid(A,B,DT,age_scale);          %make temporal grid points
CalCurves = read_curves(X,A,B,age_scale); %extract calibration curves for 14C ag[T,D0,S,phi,sigma2] = initialization(Y,F_14C,F_err,A,B,dt,year_scale,cal_dataset); %initialize parameters
%% run MCMC simulations
M = length(X.age);      %number of ages to be modeled 
N = length(G);          %number of sea-level data points to be inferred 
T_mcmcsamples = zeros(nsamples,M,nchains);
S_mcmcsamples = zeros(nsamples,N,nchains);
phi_mcmcsamples = zeros(nsamples,N,nchains);
sigma2_samples = zeros(nsamples,nchains);
kappa2_samples = zeros(nsamples,nchains);
disp('Begin Bayesian sea-level modeling...');
parpool('local');
parfor i = 1:nchains
    disp(['Generating chain ' num2str(i) '...']);
    T = initial_T(X,CalCurves,age_scale,G);
    [T_mcmcsamples(:,:,i),S_mcmcsamples(:,:,i),phi_mcmcsamples(:,:,i),sigma2_samples(:,i),kappa2_samples(:,i)] = mcmc(T,X,Y,age_scale,CalCurves,G,nsamples,burn_in,thin);
end
delete(gcp('nocreate'));
disp('Bayesian sea-level modeling completed successfully!');
%% monitor convergence
disp('Monitor the convergence of chains...')
convergence(T_mcmcsamples,'T');
convergence(S_mcmcsamples,'S');
convergence(phi_mcmcsamples,'phi');
%% mix the chains and save the results 
disp('Post processing and plotting the results...');
T_samples = round(mean(T_mcmcsamples,3)); 
S_samples = mean(S_mcmcsamples,3);
phi_samples = mean(phi_mcmcsamples,3);
sigma2_samples = mean(sigma2_samples,2);
kappa2_samples = mean(kappa2_samples,2);
save_mcmc(T_samples,S_samples,phi_samples,X,G,age_scale);
%% generate descriptive statistics and save the results
[CAL_AGE,T_pdfs,RSL,Rate] = post_process(T_samples,S_samples,phi_samples,X,Y,age_scale,A,B,G,delta,alpha);
%% plot results
figure(1)
%plot SLIPs against unmodeled ages
plot_slips0(X,Y,CalCurves,age_scale,A,B);
%
figure(2)
%plot modeld sea-level curve
plot_rsl(RSL,age_scale,Y,A,B);
hold on
%plot SLIPs against modeled ages
plot_slips1(T_pdfs,CAL_AGE,Y,A,B,age_scale);
%
figure(3) 
plot_rate(Rate,age_scale,Y,A,B); % plot rate of rsl changes
%% clean up the workspace
clear M N DT i alpha sigma2 nsamples nchains burn_in thin delta ...
    T S phi phi_mcmcsamples S_mcmcsamples T_mcmcsamples;