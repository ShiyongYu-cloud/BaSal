function X = data_process(X,age_scale)
%% function for converting 14C ages to F14C and non-14C ages to the designated
%age scale
%INPUT
%X: structure containing original info about the ages of SLIPs
%age_scale: scale of years (BC/AD, BP, B2K)
%OUTPUT
%X: structure containing augumented info about the ages of SLIPs
%%
N = length(X.age);
%correct for reservoir effect, if any, and convert to F14C space
for i = 1:N
    if strcmpi(X.age_type(i),'14C') == 1        %14C ages
        [F_14C,F_err] = radiocarbon2f(X.age(i),X.age_err(i),X.r_age(i),X.r_err(i));
        X.F14C_val(i) = F_14C;
        X.F14C_err(i) = F_err;
    elseif strcmpi(X.age_type(i),'14C') == 0    %non-14C ages
        X.F14C_val(i) = NaN;
        X.F14C_err(i) = NaN;
    end 
end    
%convert non-14C ages to the common age scale
for i = 1:N
    if strcmpi(X.age_type(i),'14C') == 1        %14C ages
        X.age_common(i) = NaN;
    elseif strcmpi(X.age_type(i),'14C') == 0    %non-14C ages  
        if strcmpi(age_scale,'BC/AD') == 1 
            X.age_common(i) = X.year_dated(i) - X.age(i);
        elseif strcmpi(age_scale,'BP') == 1 
            X.age_common(i) = X.age(i) - (X.year_dated(i) - 1950);
        elseif strcmpi(age_scale,'B2K') == 1 
            X.age_common(i) = X.age(i) - (X.year_dated(i) - 2000);     
        end
    end
end    
return;
%%
function [F_14C,F_err] = radiocarbon2f(age,age_err,r_age,r_err)
%% function for correcting reservoir effect and converting radiocarbon ages to F14C
%INPUT
%age: Laboratory radiocarbon ages
%age_err: 1 sigma standard deviation of radiocarbon ages
%r_age: local reservoir age or reservoir age offset
%r_err: 1 sigma standard deviation of local reservoir ages or age offsets
%OUTPUT
%F_14C: 14C age in the F14C space
%F_err: one standard deviaton of F14C
%%
%% correct for the reservoir effect
age = age - r_age;
age_err = sqrt(age_err.^2 + r_err.^2);
%% convert to F14C space
F_14C = exp(age/-8033);
F_err = F_14C.*age_err/8033;
return;