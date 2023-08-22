function [CAL_AGE,T_pdfs,RSL,Rate] = post_process(T_samples,S_samples,phi_samples,X,Y,age_scale,A,B,G,delta,alpha)
%% function for generating descriptive statistics of T, S, and phi from MCMC runs
%INPUT
% T_samples: calendar ages generated from MCMC run
% S_samplesL RSL generated from MCMC run
% phi_samples: rate of RSL changes generated from MCMC run
% X: structure containing the age info of SLIPs
% Y: structure containing the altitudinal info of SLIPs
% age_scale: scale of the modeled calendar ages (BC/AD, BP, B2K)
% A: early boundary of the study period
% B: late boundary of the study period
% G: temporal grid points
% delta: time resolution of the pdf or the nearest years to be rounded
% alpha: significance level
%OUTPUT
% CAL_AGE: structure containing modeled calendar ages
% T_pdfs: empirical pdf of the modeled calendar ages 
% RSL: structure containing descriptive statistics of modeled relative sea level
% Rate: structure containing descriptive statistics of modeled rate of RSL changes
% rsl.txt: data file of modeled relative sea level
% rate.txt: data file of modeled rate of rsl changes
%%
%% build up the empirical pdf of the modeled calendar ages and save the results to files
T_pdfs = mcmc2pdf(T_samples,X,A,B,delta,age_scale);
%% calculate the highest posterior density regions of the modeled calendar ages
CAL_AGE = pdf2hpd(T_pdfs,delta,X,age_scale); 
%% generate descriptive statistics of the modeled RSL
RSL = struct('posterior_mean',[],'cal_age',[],'CI_upper',[],'CI_lower',[]);
Rate = struct('posterior_mean',[],'cal_age',[],'CI_upper',[],'CI_lower',[]);
%% RSL
RSL.cal_age = G;
RSL.posterior_mean = mean(S_samples)';  
M = size(S_samples,2);
S_CI = zeros(M,2);
for j = 1:M
    S_CI(j,:) = prctile(S_samples(:,j),abs([0,100]-(100-100*(1-alpha/2))));
end
RSL.CI_lower = S_CI(:,1);
RSL.CI_upper = S_CI(:,2);
% save results
unit = char(unique(Y.unit));
header1 = cell(1,4);
header1{1,1} = strcat('Cal_age','(',age_scale,')');
header1{1,2} = strcat('Mean SL','(',unit,')');
header1{1,3} = strcat(num2str(100*(1-alpha)),'%','CI_lower','(',unit,')');
header1{1,4} = strcat(num2str(100*(1-alpha)),'%','CI_upper','(',unit,')');
fmt = repmat('%20g ',1,3); 
fmt = ['%20d ' fmt];
fid1 = fopen('rsl.txt','wt');
fprintf(fid1,'%20s ',header1{1:end});
fprintf(fid1,'\n');
for k = 1:M
    fprintf(fid1,fmt,[RSL.cal_age(k,:) RSL.posterior_mean(k,:) RSL.CI_lower(k,:) RSL.CI_upper(k,:)]);
    fprintf(fid1,'\n');
end    
fclose(fid1);
%% rate of RSL changes
Rate.cal_age = G;
Rate.posterior_mean = mean(phi_samples)';
N = size(phi_samples,2);
B_CI = zeros(N,2);
for j = 1:N
    B_CI(j,:) = prctile(phi_samples(:,j),abs([0,100]-(100-100*(1-alpha/2))));
end
Rate.CI_lower = B_CI(:,1);
Rate.CI_upper = B_CI(:,2);
% save results
%convert to "mm" for rate of rsl changes  
if strcmpi(unit,'m') == 1
    unit = 'mm';
    Rate.posterior_mean = 1000*Rate.posterior_mean;
    Rate.CI_lower = 1000*Rate.CI_lower;
    Rate.CI_upper = 1000*Rate.CI_upper;
elseif strcmpi(unit,'cm') == 1
    unit = 'mm';
    Rate.posterior_mean = 10*Rate.posterior_mean;
    Rate.CI_lower = 10*Rate.CI_lower;
    Rate.CI_upper = 10*Rate.CI_upper;
end    
% 
header2 = cell(1,4);
header2{1,1} = strcat('Cal_age','(',age_scale,')');
header2{1,2} = strcat('Mean rate','(',unit, '/yr',')');
header2{1,3} = strcat(num2str(100*(1-alpha)),'%','CI_lower','(',unit, '/yr',')');
header2{1,4} = strcat(num2str(100*(1-alpha)),'%','CI_upper','(',unit, '/yr',')');
fmt = repmat('%20g ',1,3); 
fmt = ['%20d ' fmt];
fid2 = fopen('rate.txt','wt');
fprintf(fid2,'%20s ',header2{1:end});
fprintf(fid2,'\n');
for k = 1:N
    fprintf(fid2,fmt,[Rate.cal_age(k,:) Rate.posterior_mean(k,:) Rate.CI_lower(k,:) Rate.CI_upper(k,:)]);
    fprintf(fid2,'\n');
end    
fclose(fid2);
end