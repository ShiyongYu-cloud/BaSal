function epdfs = mcmc2pdf(T_samples,X,A,B,delta,age_scale)
%% function for estimating the empirical pdf of the modeled ages and save the results 
%INPUT
% T_samples: calendar ages of the SLIPs generated from MCMC run
% X: structure containing information about ages of the SLIPs
% A: early boundary of the age sequence
% B: late boundary of the age sequence
% delta: time resolution of the empirical pdf
% age_scale: scale of year to be reported (BC/AD, BP, B2K)
%OUTPUT
% epdfs: array containing the empirical pdf of modeled calendar ages 
%%
%% preallocate size
if strcmpi(age_scale,'BC/AD') == 1
   t = A:1:B;
else
   t = B:1:A;
end
%% resample at delta-yr intervals
t = t(mod(t,delta) == 0);
M = length(t);
N = size(T_samples,2);
epdfs = zeros(M,N+1);
hpd68_2 = zeros(M,N+1);
hpd95_4 = zeros(M,N+1);
%% get the pdfs and the one- and two-sigma posterior probability of the modeled ages
for i = 1:N
    [cal_age,prob,hp_68p2,hp_95p4] = epdf(T_samples(:,i),A,B,delta,age_scale);
    epdfs(:,i+1) = prob;
    hpd68_2(:,i+1) = hp_68p2;
    hpd95_4(:,i+1) = hp_95p4;
end
epdfs(:,1) = cal_age;
hpd68_2(:,1) = cal_age; 
hpd95_4(:,1) = cal_age;
%% save results to a file
outfiles = {'T_pdfs.txt','T_hp68_2.txt','T_hp95_4.txt'};
results = [epdfs hpd68_2 hpd95_4];
for i = 1:length(outfiles)
    save_pdf(outfiles{i},X,age_scale,results(:,(N+1)*i-N:(N+1)*i));
end
return;
%%
function [cal_age,prob,hp_68p2,hp_95p4] = epdf(Z,A,B,delta,age_scale)
%% function for estimating the empirical pdf and the highest posterior probability of the calendar ages
%INPUT
% Z: vector containing random numbers obtained from MCMC run 
% A: early boundary of the study period
% B: late boundary of the study period
% delta: time resolution of the empirical pdf in years
% age_scale: scale of year(BC/AD, BP, B2K)
%OUTPUT
% cal_age: vector containing age points  
% prob: vector containing the probability of the age points
% hp_68p2: 68.2% posterior probability of the age points
% hp_95p4: 95.4% posterior probability of the age points
%% building up the empirical pdf
%generate points at which the pdf will be estimated
if strcmpi(age_scale,'BC/AD') == 1
   cal_age = A:1:B;
else
   cal_age = B:1:A;
end
%resample at delta-yr interval
cal_age = cal_age(mod(cal_age,delta) == 0);
cal_age = cal_age';
M = length(cal_age);
N = length(Z);
prob = zeros(M,1);
hp_68p2 = zeros(M,1);
hp_95p4 = zeros(M,1);
F = zeros(M,1);
h = 0.000000001;        % Small value closer to zero for evaluating
                        % numerical differentiation.
% Estimating CDF by its definition
for i = 1:M
    p = 0;              % True Probability
    q = 0;              % False Probability
    for j = 1:N
        if Z(j) <= cal_age(i)   % Definition of CDF
            p = p + 1;
        else
            q = q + 1;
        end
    end
    F(i) = p/(p + q);   % Calulating Probability
end
% Estimating PDF by differentiating the CDF
for k = 1:M
    fxph = interp1(cal_age,F,cal_age(k) + h,'spline');  % Interpolating value of F(x+h)
    fxmh = interp1(cal_age,F,cal_age(k) - h,'spline');  % Interpolating value of F(x-h)
    prob(k) = (fxph - fxmh)/(2*h); 
    if prob(k) < 0
       prob(k) = 0;
    end        
end                                         
prob = smooth(prob);  
% Normalizing to 1
prob = prob./sum(prob);
%% calculate the highest posterior probabilities
calprob = [cal_age prob];
hpd = calprob(:,1:2);
hpd = sortrows(hpd, 2);
hpd(:,3) = cumsum(hpd(:,2));
%One sigma posterior probability
hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)), :);
hpd68_2 = sortrows(hpd68_2,1);
id0 = ismember(cal_age,hpd68_2(:,1)); 
hp_68p2(id0) = hpd68_2(:,2);
%two sigma posterior probability
hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)), :);
hpd95_4 = sortrows(hpd95_4,1);
id1 = ismember(cal_age,hpd95_4(:,1)); 
hp_95p4(id1) = hpd95_4(:,2);
return;
%%
function [] = save_pdf(file_name,X,age_scale,pdfs)
%% function for saving the empirical pdfs to a file
%INPUT
% file_name: file name for data to be saved
% X: structure containing information about ages of the SLIPs
% age_scale: scale of year (BC/AD, BP, B2K)
% pdfs: probability of the calendar ages
%% creat a header
N = length(X.age); 
header = cell(1,N+1); 
header{1} = strcat('Age (',age_scale,')');
for i = 1:N
    header{1,1+i} = X.sample_ID{i};
end
%% write header and data
fmt = repmat('%25g ',1,N); 
fmt = ['%25d ' fmt];
fid = fopen(file_name,'wt');
fprintf(fid,'%25s ',header{1:end});
fprintf(fid,'\n');
for k = 1:size(pdfs,1)
    fprintf(fid,fmt,pdfs(k,:));
    fprintf(fid,'\n');
end    
fclose(fid);
return;