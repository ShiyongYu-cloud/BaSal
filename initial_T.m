function T = initial_T(X,CalCurves,age_scale,G)
%% function for initializing model parameters
%INPUT
%X: structure containing augumented info about the ages of SLIPs
%CalCurve: sturcture containing calibration curves for 14C ages in BP scale
%age_scale: scale of year to be reported (e.g. BC/AD, BP, or B2K)
%G: temporal grid points
%OUTPUT
%T: initial values of the calendar age sequence
%%
%% obtain the initial values of the calendar ages 
A = G(1);
B = G(end);
cal_age = gen_age(X,CalCurves,age_scale,A,B);  
if strcmpi(age_scale,'BC/AD') == 1 
    true = all(cal_age > A) && all(cal_age < B) && length(unique(cal_age)) == length(cal_age); % unique ages 
    while(~true)
         cal_age = gen_age(X,CalCurves,age_scale,A,B);   
         true = all(cal_age > A) && all(cal_age < B) && length(unique(cal_age)) == length(cal_age);
    end
    T = sort(cal_age,'ascend'); % make all(diff(cal_age))>0
else
    true = all(cal_age < A) && all(cal_age > B) && length(unique(cal_age)) == length(cal_age);% unique ages 
    while(~true)
        cal_age = gen_age(X,CalCurves,age_scale,A,B);   
        true = all(cal_age < A) && all(cal_age > B) && length(unique(cal_age)) == length(cal_age);
    end
    T = sort(cal_age,'descend'); % make all(diff(cal_age))<0
end
return;
%%
function cal_age = gen_age(X,CalCurves,age_scale,A,B)
%% function for generating a calendar age sequence 
%INPUT
%X: structure containing information about ages of the SLIPs
%CalCurve: structure containing the calibration curves
%age_scale: scale of year (e.g. BC/AD, BP, B2K)
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%OUTPUT
%cal_age: a vector containing the generated calendar ages
%%
N = length(X.age);
cal_age = zeros(N,1);
dt = 1; % interpolate the curve to annual resolution
%
for i = 1:N
    if strcmpi(X.age_type(i),'14C') == 0       %non-14C age
        %draw a random number from the truncated normal distribution
        norm_dist = makedist('Normal',X.age_common(i),X.age_err(i));
        if strcmpi(age_scale,'BC/AD') == 1 
            truncated_norm = truncate(norm_dist,A,B);
        else  
            truncated_norm = truncate(norm_dist,B,A);
        end    
        x = random(truncated_norm,1,1);
        cal_age(i) = round(x);
    elseif strcmpi(X.age_type(i),'14C') == 1   %14C age 
        % retrieve the calibration dataset
        curve_cal_age = CalCurves.(char(X.cal_curve(i)))(:,1);
        curve_F14_val = CalCurves.(char(X.cal_curve(i)))(:,2);
        curve_F14_err = CalCurves.(char(X.cal_curve(i)))(:,3);
        hi_curve_cal_age = curve_cal_age(1):dt:curve_cal_age(end);
        hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
        hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
        % calculate the highest posterior probabilities
        a = (X.F14C_val(i) - hi_curve_F14_val).^2;
        b = 2*(X.F14C_err(i)^2 + hi_curve_F14_err.^2);
        c = sqrt(X.F14C_err(i)^2 + hi_curve_F14_err.^2);
        cal_prob = exp(-a./b)./c;
        % normalize the probabilities to 1
        cal_prob = cal_prob/sum(cal_prob);
        % draw a calendar age from the pdf
        x = durand(hi_curve_cal_age,cal_prob,1,1);
        cal_age(i) = round(x);
        % convert to the designated age scale
        if strcmpi(age_scale,'BC/AD') == 1
            cal_age(i) = 1950 - cal_age(i);
        elseif strcmpi(age_scale,'B2K') == 1
            cal_age(i) = 50 + cal_age(i);
        end
    end
end
return;
%%
function x = durand(values,prob,M,N)
%% function for drawing M*N random numbers from a list at given probability
%INPUT:
%values: a list of numbers to be drawn
%prob: corresponding probability value of the numbers
%M: number of rows of random numbers 
%N: number of columns of random numbers
%OUTPUT:
%x: M*N random numbers  
%%
values = values(:);
prob = prob(:);
if sum(prob)~=1
   prob = prob/sum(prob);
end
L = length(prob);
K = M*N;
psup = cumsum(prob);
pinf = [0; psup(1:end-1)];
Pinf = kron(ones(1,K),pinf(:));
Psup = kron(ones(1,K),psup(:));
u = rand(1,K);
U = kron(ones(L,1),u);
C = (U>Pinf) & (U<Psup);
V = kron(values(:),ones(1,K));
X = V.*C;
x = sum(X);
x = reshape(x,M,N);
return;