function [] = plot_rsl(RSL,age_scale,Y,A,B)
%% function for plotting modeled RSL curve
%INPUT
% RSL: structure containing descriptive statistics of the modeled RSL
% age_scale: scale of year to be reported (e.g. BC/AD, BP, B2K)
% Y: structure containing the altitudinal information of the SLIPs
% A: early boundary of the study period
% B: late boundary of the study period
%%
%% plot the modeled sea-level curve
cal_age = RSL.cal_age(1:end);
post_mean = RSL.posterior_mean(1:end);
CI_lower = RSL.CI_lower(1:end);
CI_upper = RSL.CI_upper(1:end);
XX = [cal_age; flipud(cal_age); cal_age(1)];
YY = [CI_lower; flipud(CI_upper); CI_lower(1)];
%XX = [RSL.cal_age; flipud(RSL.cal_age); RSL.cal_age(1)];
%YY = [RSL.CI_lower; flipud(RSL.CI_upper); RSL.CI_lower(1)];
%h1 = fill(XX,YY,[0.8 0.8 0.8],'LineStyle','none');
%h1 = fill(XX,YY,[91, 207, 244]/255,'LineStyle','none');
h = fill(XX,YY,'g','LineStyle','none');
set(h,'facealpha',0.5);
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
set(gca, 'TickDir', 'out');
if strcmpi(age_scale,'BP') == 1 
   set(gca,'XDir','reverse');
   xlabel('Calendar age (BP)');
elseif strcmpi(age_scale,'B2K') == 1 
   set(gca,'XDir','reverse');
   xlabel('Calendar age (B2K)');
else
   xlabel('Calendar age (BC/AD)');    
end
unit = char(unique(Y.unit)); 
ylabel(strcat('RSL (', unit,')'));
grid on;
%legend('95% confidence interval','Mean sea level','location','southeast'); 
end