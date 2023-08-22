function [X,Y] = read_data(data_file)
%% function for reading the SLIPs from a data file
%INPUT
%data_file: data file of the SLIPs
%OUTPUT
%X: structure containing info about the ages of SLIPs
%Y: structure containing info about elevation of SLIPs
%%
%% reading data
h = fopen(data_file);
DATA = textscan(h,'%s %f %f %s %f %f %s %f %f %s %f' ,'headerlines',1,'delimiter','\t');
fclose(h);
sample_ID = DATA{1};
rsl = DATA{2};
rsl_err = DATA{3};
unit = DATA{4};
age = DATA{5};
age_err = DATA{6};
age_type = DATA{7};
r_age = DATA{8};
r_err = DATA{9};
cal_curve = DATA{10};
year_dated = DATA{11};
%convert 'NaA' to 0
r_age(isnan(r_age)) = 0;
r_err(isnan(r_err)) = 0;
%% wrapping data into structures
X = struct('sample_ID',[],'age',[],'age_err',[],'age_type',[],'r_age',[],'r_err',[],'cal_curve',[],'year_dated',[]);
X.sample_ID = sample_ID;
X.age = age;
X.age_err = age_err;
X.age_type = age_type;
X.r_age = r_age;
X.r_err = r_err;
X.cal_curve = cal_curve;
X.year_dated = year_dated;
%
Y = struct('rsl',[],'rsl_err',[],'unit',[]);
Y.rsl = rsl;
Y.rsl_err = rsl_err;
Y.unit = unit;
end