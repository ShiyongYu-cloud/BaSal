%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BaSaL: A MATLAB package for inferring the pattern and rate of past     % 
%              RSL changes                                               %
% Copyright: Shiyong Yu                                                  %  
% E-mail: shiyong.yu@gmail.com                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clc;
clear;
%% synthesize data
tau = 15000:-100:0; % calendar ages
y = (tau-12500).*(tau-7500).*tau.^2*1e-14; % mean sea level
% randomly draw 20 calendar ages from [300 14800] 
x = 300 + random('unid',14500,20,1); 
x = sort(x,'descend'); 
x = round(x/100)*100;
x = unique(x,'stable'); 
% randomly pick 10 calendar 14C ages 
r_id = sort(randperm(numel(x),10),'ascend'); 
r_cal = x(r_id);  
% convert calendar 14C ages to conventional 14C ages
[c14_age,c14_err] = cal2radiocarbon(r_cal,'intcal20');
% wrap up the data
X = struct('sample_ID',[],'age',[],'age_err',[],'age_type',[],'r_age',[],'r_err',[],'cal_curve',[],'year_dated',[]);
Y = struct('rsl',[],'rsl_err',[],'unit',[]);
X.sample_ID = cell(numel(x),1);
X.age_type = cell(numel(x),1);
X.cal_curve = cell(numel(x),1);
for i = 1:numel(x)
    X.sample_ID{i} = strcat('Sample_',num2str(i));
    X.age_type{i} = 'OSL';
    X.cal_curve{i} = 'None';
end
X.age = x;
X.age(r_id) = round(c14_age/10)*10;
X.age_err = round(random('gam',50,2,numel(x),1)/10)*10;
X.age_err(r_id) = round(c14_err/10)*10;
X.age_type(r_id) = cellstr('14C');
X.cal_curve(r_id) = cellstr('intcal20');
X.r_age = zeros(numel(x),1);
X.r_err = zeros(numel(x),1);
X.year_dated = 2000*ones(numel(x),1);
y_mean = interp1(tau,y,x);
Y.rsl_err = random('gam',2,1/2,numel(x),1);
Y.unit = cell(numel(x),1);
for i = 1:numel(x)
    pd = makedist('normal',y_mean(i),Y.rsl_err(i));
    tn = truncate(pd,y_mean(i)-Y.rsl_err(i),y_mean(i)+Y.rsl_err(i));
    Y.rsl(i) = random(tn,1,1);
    Y.unit{i} = 'm';
end    
%% set up model parameters
age_scale = 'BP';                 % specify age scale (BC/AD,BP,B2K)
DT = 25;                           % specify the size of time grid
delta = 10;        % specify the nearest years the modeled calendar ages to be rounded  
alpha = 0.05;     % specify the significance level for estimating confidence interval
nsamples = 500;   % specify the number of samples to be kept for each chain
nchains = 3;      % specify the number of chains to be run for convergence diagnosis
burn_in = 1000;   % specify the burn-in period of the Markov chains
thin = 20;        % specify the thinning interval of the Markov chains
%% read and preprocess data (convert to F14C for 14C ages and to the specified age scale for non-14C ages)
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
subplot(3,1,1);
plot_slips0(X,Y,CalCurves,age_scale,A,B);
yticks([-10 0 10 20 30 40]);
ylim([-10 40]);
%
%figure(2)
%plot modeld sea-level curve
subplot(3,1,2);
plot_rsl(RSL,age_scale,Y,A,B);
hold on
%plot SLIPs against modeled ages
plot_slips1(T_pdfs,CAL_AGE,Y,A,B,age_scale);
yticks([-10 0 10 20 30 40]);
ylim([-10 40]);
%
%figure(3) 
subplot(3,1,3);
plot_rate(Rate,age_scale,Y,A,B); % plot rate of rsl changes
ylim([-20 10]);
%% clean up the workspace
clear M N DT i alpha sigma2 nsamples nchains burn_in thin delta ...
    T S phi phi_mcmcsamples S_mcmcsamples T_mcmcsamples;