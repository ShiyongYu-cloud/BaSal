function CalCurves = read_curves(X,A,B,age_scale)
%% function for reading calibration curves associated with the 14C ages 
%INPUT:
%X: structure containing information about ages of the SLIPs
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%age_scale: scale of year (e.g. BC/AD, BP, B2K)
%OUTPUT
%CalCurve: sturcture containing the calibration curves for 14C ages in BP
%scale
%%
id = strcmpi(X.age_type,'14C') == 1;  %14C ages
curves = X.cal_curve(id);
curves = unique(curves); % extract unique calibration curves
CalCurves = struct();    % put the caruves into a structture 
%% read and extract the calibration curves and put them into a structure
for i = 1:size(curves,1)
     CalCurves.(curves{i}) = extract_curve(curves{i},A,B,age_scale);
end
return;
%%
function cal_dataset = extract_curve(cal_curve,A,B,age_scale)
%% function for reading the calibration curve and converting it to the F14C space
%INPUT
%cal_curve: name of the calibration curve
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%age_scale: scale of year (e.g. BC/AD, BP, B2K)
%OUTPUT
%cal_dataset: array containing the calibration curve bounded by A and B
%% open the calibration curve data
headerlines = 11;
h = fopen([cal_curve,'.14c']);
cal_dataset = textscan(h,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(h);
curve_cal_age = flipud(cal_dataset{1});
curve_C14_age = flipud(cal_dataset{2});
curve_C14_err = flipud(cal_dataset{3});
%% convert radiocrbon ages in the calibration curve to F_14C space
curve_F14_val = exp(curve_C14_age/-8033); %convert the radiocarbon ages to the F_14C space
curve_F14_err = curve_F14_val.*curve_C14_err/8033; %convert the error to the F_14C space
%% convert the early and late boundaries to BP scale
if strcmpi(age_scale,'BC/AD') == 1
    A = 1950 - A;
    B = 1950 - B;
elseif strcmpi(age_scale,'B2K') == 1
    A = A - 50;
    B = B - 50;   
end
%% extract calibration dataset
ind = (curve_cal_age >= B & curve_cal_age <= A);
cal_dataset = [curve_cal_age(ind) curve_F14_val(ind) curve_F14_err(ind)];
return;