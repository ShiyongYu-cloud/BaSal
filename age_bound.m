function [A,B] = age_bound(X,age_scale)
%% function for estimating the early and late boundary of the age sequence
%INPUT
% X: structure containing information about ages of the SLIPs
% age_scale: scale of year to be reported (e.g. BC/AD, BP, or B2K)
%OUTPUT
% A: early boundary of the age sequence
% B: late boundary of the age sequence
%%
p_threshold = 0.0005;
% find the early bound of the age sequence
if strcmpi(X.age_type(1),'14C') == 0            %non_14C age
    A = X.age_common(1);
    pA = normpdf(A,X.age_common(1),X.age_err(1));
    while pA > p_threshold
        if strcmpi(age_scale,'BC/AD') == 1 
            A = A - 10;
        else
            A = A + 10;
        end    
        pA = normpdf(A,X.age_common(1),X.age_err(1));
    end
elseif strcmpi(X.age_type(1),'14C') == 1        %14C age
    A_F_14C = X.F14C_val(1);
    A_F_err = X.F14C_err(1);
    %find the associated calibration curve
    A_cal_curve = X.cal_curve{1};
    % calibrate the age
    [A_cal_age,A_prob] = calibration(A_F_14C,A_F_err,A_cal_curve);
    A = max(A_cal_age(A_prob > p_threshold));
    % conver to other age scales
    if strcmpi(age_scale,'BC/AD') == 1
        A = 1950 - A;
    elseif strcmpi(age_scale,'B2K') == 1
        A = 50 + A;
    end
end
%round to the nearest 100 years
if strcmpi(age_scale,'BC/AD') == 1
    A = floor(A/100)*100; 
else
    A = ceil(A/100)*100; 
end
% find the late bound of the age sequence
if strcmpi(X.age_type(end),'14C') == 0          %non-14C age
    B = X.age_common(end);
    pB = normpdf(B,X.age_common(end),X.age_err(end));
    while pB > p_threshold
        if strcmpi(age_scale,'BC/AD') == 1 
            B = B + 10;
        else
            B = B - 10;
        end    
        pB = normpdf(B,X.age_common(end),X.age_err(end));
    end
elseif strcmpi(X.age_type(end),'14C') == 1      %14C age
    B_F_14C = X.F14C_val(end);
    B_F_err = X.F14C_err(end);
    %find the associated calibration curve
    B_cal_curve = X.cal_curve{end};
    % calibrate the age
    [B_cal_age,B_prob] = calibration(B_F_14C,B_F_err,B_cal_curve);
    B = min(B_cal_age(B_prob > p_threshold));
    % conver to other age scales
    if strcmpi(age_scale,'BC/AD') == 1
       B = 1950 - B;
    elseif strcmpi(age_scale,'B2K') == 1
       B = 50 + B;
    end
end
%round to the nearest 100 years
if strcmpi(age_scale,'BC/AD') == 1
    B = ceil(B/100)*100; 
    if B > 2000
        B = 0;
    end    
else
    B = floor(B/100)*100; 
    if B < 0
        B = 0;
    end    
end
return;
%%
function [cal_age,prob] = calibration(F_14C_val,F_14C_err,cal_curve)
%% function for calibrating a 14C age using Bayesian higher posterior density analysis
%INPUT
% F_14C:    Lab 14C age in F14C space
% F_err:    Lab 14C uncertainty in F14C space
% cal_curve: name of the calibration curve
%OUTPUT
% cal_age:  calibrated age in year_scale (BC/AD, BP, B2K)
% prob: probability of the calibrated ages
%%
% Please see manuscript for more detailed information:
% Lougheed, B.C. & Obrochta, S.P. (2016). MatCal: Open Source Bayesian
% 14C Age Calibration in Matlab. Journal of Open Research Software. 4(1),
% p.e42. DOI: http://doi.org/10.5334/jors.130
%%
%% open the calibration curve data file
headerlines = 11;
h = fopen([cal_curve,'.14c']);
cal_dataset = textscan(h,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(h);
curve_cal_age = flipud(cal_dataset{1});
curve_C14_age = flipud(cal_dataset{2});
curve_C14_err = flipud(cal_dataset{3});
%% convert radiocrbon ages in the calibration curve to F14C space
curve_F14_val = exp(curve_C14_age/-8033); %ages to F14C
curve_F14_err = curve_F14_val.*curve_C14_err/8033; %errors to F14C
% interpolate the calendar age in the curve to annual resolution
dt = 1;
hi_curve_cal_age = curve_cal_age(1):dt:curve_cal_age(end);
hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
%% Calculate probability for every calendar year in the F14 space
cal_age = hi_curve_cal_age;
a = (F_14C_val - hi_curve_F14_val).^2;
b = 2*(F_14C_err^2 + hi_curve_F14_err.^2);
c = sqrt(F_14C_err^2 + hi_curve_F14_err.^2);
prob = exp(-a./b)./c;
%% normalize the probabilities to 1
prob = prob/sum(prob);
return;