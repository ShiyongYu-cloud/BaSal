function [] = plot_comparison(T_pdfs,CAL_AGE,RSL,Rate,age_scale,Y,A,B,IGP)
%% function for plotting the results compared with the EIV-IGP model
%INPUT
% T_pdfs: empirical pdf of the modeled calendar ages 
% CAL_AGE: structure containing modeled calendar ages
% RSL: structure containing the modeled RSL
% Rate: structure containing the modeled rates of RSL changes
% age_scale: scale of year to be reported (e.g. BC/AD or BP)
% Y: structure containing the altitudinal information of the SLIPs
% A: early bound of the study period
% B: late bound of the study period
% IGP: IGP results
%%
%subplot(2,1,1); 
%plot rsl curve
plot_rsl(RSL,age_scale,Y,A,B);
hold on
%plot SLIPs against modeled ages
plot_slips1(T_pdfs,CAL_AGE,Y,A,B,age_scale);
hold on
%plot EIV-IGP results
errorbar(IGP(:,1),IGP(:,2),IGP(:,3));
set(gca,'XDir','reverse');
%%
%subplot(2,1,2); 
%plot rate of rsl changes
plot_rate(Rate,age_scale,Y,A,B);
hold on
%plot EIV-IGP results
errorbar(IGP(:,1),IGP(:,4)/1000,IGP(:,5)/1000);
set(gca,'XDir','reverse');
end