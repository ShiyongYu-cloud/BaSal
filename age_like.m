function p = age_like(X,T,A,B,age_scale,CalCurves)
%% function for calculating the likelihood given calendar ages
%INPUT
%X: structure containing information about ages of the SLIPs
%T: modeled calendar ages
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%age_scale: scale of year (BC/AD, BP, B2K)
%CalCurves: structure containing the calibration curves for 14C ages
%OUTPUT
%p: likelihood value in log scale
%%
%% deal with non-14C ages
%extract non_14C ages
id_cal = strcmpi(X.age_type,'14C') == 0;
model_cal_age = T(id_cal);
if isempty(model_cal_age)
    p_cal = log(1);
else    
    if strcmpi(age_scale,'BC/AD') == 1
        if any(model_cal_age > A) || any(model_cal_age < B) 
           p_cal = log(1);
        else
           p_cal = log(0);
        end   
    else
        if any(model_cal_age < A) || any(model_cal_age > B) 
           p_cal = log(1);
        else
           p_cal = log(0); 
        end
    end
    lab_cal_age = X.age_common(id_cal)';
    lab_cal_err = X.age_err(id_cal);
    a = (lab_cal_age - model_cal_age).^2;
    b = 2*lab_cal_err.^2;
    c = sqrt(2*pi*lab_cal_err.^2);
    p_cal = p_cal + sum(-log(c))+sum(-a./b);
end
%% deal with 14C ages
%extract 14C ages
id_rad = strcmpi(X.age_type,'14C') == 1;
model_rad_age = T(id_rad);
if isempty(model_rad_age)
    p_rad = log(1);
else    
    if strcmpi(age_scale,'BC/AD') == 1
        if any(model_rad_age > A) || any(model_rad_age < B) 
           p_rad = log(1);
        else
           p_rad = log(0);
        end   
    else
        if any(model_rad_age < A) || any(model_rad_age > B) 
           p_rad = log(1);
        else
           p_rad = log(0); 
        end
    end
    lab_F14C_val = X.F14C_val(id_rad)';
    lab_F14C_err = X.F14C_err(id_rad)';
    cal_curves = X.cal_curve(id_rad);
    N = length(model_rad_age);
    %find the F14C and errors in the calibration curve 
    [model_F14C_val, model_F14C_err] = calendar2f(model_rad_age,cal_curves,CalCurves,age_scale);
    a = (lab_F14C_val - model_F14C_val).^2;
    b = 2*(lab_F14C_err.^2 + model_F14C_err.^2);
    c = gamma(N/2);
    %c = sqrt(2*pi*(lab_F14C_err.^2 + model_F14C_err.^2));
    %p_rad = p_rad + sum(-log(c))+sum(-a./b);
    p_rad = p_rad + log(c)-(N/2)*log(sum(a./b));
end
p = p_cal + p_rad;
return;
%%
function [F14_val, F14_err] = calendar2f(cal_age,cal_curves,CalCurves,age_scale)
%% function for mapping a calendar age to radiocarbon age and converting it to F14C
%INPUT
%cal_age: calendar ages
%cal_curves: calibration curves associated with the 14C ages
%CalCurves: structure containing the extracted calibration curves
%age_scale: scale of years to be reported (BC/AD, BP, B2K)
%OUTPUT
%F14_val: F14C values of the calendar ages in the calibration curve
%F14_err: F14C errors of the calendar ages in the calibration curve  
%% convert to BP scale
if strcmpi(age_scale,'BC/AD') == 1
   cal_age = 1950 - cal_age; 
elseif strcmpi(age_scale,'B2K') == 1
   cal_age = cal_age - 50; 
end
%% find the corresponding F14C values in the calibration curve
N = length(cal_age);
F14_val = zeros(N,1);
F14_err = zeros(N,1);
for i = 1:N
    % retrieve the calibration dataset
    curve_cal_age = CalCurves.(cal_curves{i})(:,1);
    curve_F14_val = CalCurves.(cal_curves{i})(:,2);
    curve_F14_err = CalCurves.(cal_curves{i})(:,3);
    F14_val(i) = interp1(curve_cal_age, curve_F14_val, cal_age(i));
    F14_err(i) = interp1(curve_cal_age, curve_F14_err, cal_age(i));
end
return;