function G = make_grid(A,B,DT,age_scale)
%% function for generating temporal grid points
%INPUT
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%DT: size of the temporal grids
%age_scale: scale of year to be reported (e.g. BC/AD, BP, or B2K)
%OUTPUT
%G: tempora lgrid points
%%
if strcmpi(age_scale,'BC/AD') == 1
   G = A:DT:B;
else
   G = A:-DT:B;
end
G = G';
end