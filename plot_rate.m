function [] = plot_rate(Rate,age_scale,Y,A,B)
%% function for plotting the modeled rate of RSL changes
%INPUT
% Rate: structure containing descriptive statistics of the modeled rates of RSL changes
% age_scale: scale of year to be reported (e.g. BC/AD, BP, B2K)
% Y: structure containing the altitudinal information of the SLIPs
% A: early boundary of the study period
% B: late boundary of the study period
%%
%% plot modeled rate of sea-level changes
cal_age = Rate.cal_age(1:end);
post_mean = Rate.posterior_mean(1:end);
CI_lower = Rate.CI_lower(1:end);
CI_upper = Rate.CI_upper(1:end);
XX=[cal_age; flipud(cal_age)];
YY = [CI_lower; flipud(CI_upper)];
%h = fill(XX,YY,[0.8 0.8 0.8]);
h = fill(XX,YY,[91, 207, 244]/255,'LineStyle','none');
set(h,'facealpha',0.8);
hold on
plot(cal_age,post_mean,'b-','LineWidth',1);
%% set axes
set(gca,'XMinorTick','on','YMinorTick','on');
if strcmpi(age_scale,'BC/AD') == 1
    A = floor(A/500)*500;
    B = ceil(B/500)*500;
    xlim([A B]);
else
    A = ceil(A/500)*500;
    B = floor(B/500)*500;
    xlim([B A]);
end    
if strcmpi(age_scale,'BP') == 1 
   set(gca,'XDir','reverse');
   xlabel('Calendar age (BP)');
elseif strcmpi(age_scale,'B2K') == 1 
   set(gca,'XDir','reverse');
   xlabel('Calendar age (B2K)');
else
   xlabel('Calendar age (BC/AD)');  
end
%convert to "mm"
unit = char(unique(Y.unit));
if strcmpi(unit,'mm') == 0
    unit = 'mm';
end    
ylabel(strcat('Rate of RSL changes (', unit,'/yr)'));
set(gca, 'TickDir', 'out');
grid on;
legend('95 % confidence interval','Posterior mean','location', 'northeast'); 
end