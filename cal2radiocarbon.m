function [r_age,r_err] = cal2radiocarbon(cal_age,cal_curve)
%% function for converting calendar ages to conventional 14C ages
%INPUT
%cal_age: calendar ages
%cal_curve: name of the calibration curve
%OUTPUT
%r_age: conventional 14C ages
%r_err: 14C age error
%% open the calibration curve data
headerlines = 11;
h = fopen([cal_curve,'.14c']);
cal_dataset = textscan(h,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(h);
curve_cal_age = flipud(cal_dataset{1});
curve_C14_age = flipud(cal_dataset{2});
curve_C14_err = flipud(cal_dataset{3});
%% convert calendar ages to conventional 14C ages
r_age = interp1(curve_cal_age,curve_C14_age,cal_age); %convert the radiocarbon ages to the F_14C space
r_err = interp1(curve_cal_age,curve_C14_err,cal_age);
end