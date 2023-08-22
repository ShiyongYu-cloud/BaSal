function [] = plot_slips0(X,Y,CalCurves,age_scale,A,B)
%% function for plotting SLIPs against unmodeled ages
%INPUT
%X: structure containing augumented info about the ages of SLIPs
%Y: structure containing altitudinal info about the SLIPs
%CalCurve: sturcture containing calibration curves for 14C ages in BP scale
%age_scale: scale of year to be reported (e.g. BC/AD, BP, or B2K)
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%% build up pdfs of the calendar ages
cal_prob = gen_prob(X,CalCurves,age_scale,A,B);
%% calculate the highest posterior density regions of the calendar ages
delta = 1;
CAL_AGE = pdf2hpd(cal_prob,delta,X,age_scale);
%% plot the SLIPs
cal_age = cal_prob(:,1); 
prob = cal_prob(:,2:end);
[M,N] = size(prob);
prob_u = zeros(M,N);
prob_l = zeros(M,N);
for i = 1:N
    factor = Y.rsl_err(i)/max(prob(:,i));
    prob_u(:,i) = Y.rsl(i) + prob(:,i)*factor; % blow up to the upper limit of error
    prob_l(:,i) = Y.rsl(i) - prob(:,i)*factor; % blow up to the lower limit of error
    ind1 = (prob_u(:,i) > Y.rsl(i))+0.0001 & (prob_l(:,i) < Y.rsl(i)-0.0001) ;
    early = max(cal_age(ind1));
    late = min(cal_age(ind1));
    ind2 = (cal_age <= early) & (cal_age >= late);
    calage = cal_age(ind2);
    CALAGE = [calage(1); calage; calage(end); calage(1)];
    PROB_u = [Y.rsl(i); prob_u(ind2,i); Y.rsl(i); Y.rsl(i)];
    PROB_l = [Y.rsl(i); prob_l(ind2,i); Y.rsl(i); Y.rsl(i)];
    fill(CALAGE,PROB_u,[0.301 0.745 0.933],'EdgeColor','none');  % plot pdfs
    hold on
    fill(CALAGE,PROB_l,[0.301 0.745 0.933],'EdgeColor','none');  % plot pdfs
    hold on
    p95_4 = CAL_AGE(i).P95_4_Credible_intervals;    % plot the 95.4% pdf
    K = size(p95_4,1);
    for j = 1:K
        id = (cal_age >= p95_4(j,2)) & (cal_age <= p95_4(j,1));
        age = cal_age(id);
        AGE = [age(1); age; age(end); age(1)];
        pdf_u = prob_u(id,i);
        pdf_l = prob_l(id,i);
        PDF_u = [Y.rsl(i); pdf_u; Y.rsl(i); Y.rsl(i)];
        PDF_l = [Y.rsl(i); pdf_l; Y.rsl(i); Y.rsl(i)];
        fill(AGE,PDF_u,[0.929 0.694 0.125],'EdgeColor','none');
        fill(AGE,PDF_l,[0.929 0.694 0.125],'EdgeColor','none');
    end
    hold on
    p68_2 = CAL_AGE(i).P68_2_Credible_intervals;    % plot the 68.2% pdf 
    L = size(p68_2,1);
    for k = 1:L
        id = (cal_age >= p68_2(k,2)) & (cal_age <= p68_2(k,1));
        age = cal_age(id);
        AGE = [age(1); age; age(end); age(1)];
        pdf_u = prob_u(id,i);
        pdf_l = prob_l(id,i);
        PDF_u = [Y.rsl(i); pdf_u; Y.rsl(i); Y.rsl(i)];
        PDF_l = [Y.rsl(i); pdf_l; Y.rsl(i); Y.rsl(i)];
        fill(AGE,PDF_u,[0.635 0.078 0.184],'EdgeColor','none'); 
        fill(AGE,PDF_l,[0.635 0.078 0.184],'EdgeColor','none'); 
    end
    line(calage, prob_u(ind2,i),'Color','k','LineWidth',0.1);
    line(calage, prob_l(ind2,i),'Color','k','LineWidth',0.1);
end
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
unit = char(unique(Y.unit));
ylabel(strcat('RSL (', unit,')'));
set(gca, 'TickDir', 'out');
grid on;
%% make legend
x0 = (A+B)/2;
y0 = min(Y.rsl);
dy = max(Y.rsl_err)+min(Y.rsl_err);
p1.y = [y0;y0;y0+dy;y0+dy;y0];
p2.y = [y0+dy;y0+dy;y0+2*dy;y0+2*dy;y0+dy];
p3.y = [y0+2*dy;y0+2*dy;y0+3*dy;y0+3*dy;y0+2*dy];
if strcmpi(age_scale,'BC/AD') == 1
    p.x = [x0;x0+500;x0+500;x0;x0];
    fill(p.x,p1.y,[0.301 0.745 0.933]);
    text(x0+800,y0+dy/2,'Age probability','FontSize',10);
    fill(p.x,p2.y,[0.929 0.694 0.125]);
    text(x0+800,y0+1.5*dy,'95.4% HPD regions','FontSize',10);
    fill(p.x,p3.y,[0.635 0.078 0.184]);
    text(x0+800,y0+2.5*dy,'68.2% HPD regions','FontSize',10);
else
    p.x = [x0;x0-500;x0-500;x0;x0];
    fill(p.x,p1.y,[0.301 0.745 0.933]);
    text(x0-800,y0+dy/2,'Age probability','FontSize',10);
    fill(p.x,p2.y,[0.929 0.694 0.125]);
    text(x0-800,y0+1.5*dy,'95.4% HPD regions','FontSize',10);
    fill(p.x,p3.y,[0.635 0.078 0.184]);
    text(x0-800,y0+2.5*dy,'68.2% HPD regions','FontSize',10);
end
return;
%%
function cal_prob = gen_prob(X,CalCurves,age_scale,A,B)
%% function for generating the probability distribution of the calendar ages
%INPUT
%X: structure containing information about ages of the SLIPs
%CalCurve: structure containing the calibration curves
%age_scale: scale of year (e.g. BC/AD, BP, B2K)
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%OUTPUT
%cal_prob: probability distribtuion of the calendar ages
%%
N = length(X.age);
dt = 1; % interpolate the curve to annual resolution
if strcmpi(age_scale,'BC/AD') == 1
   T = A:dt:B;
else
   T = A:-dt:B;
end
T = T';
M = length(T);
cal_prob = zeros(M,N+1);
cal_prob(:,1) = T;
%
for i = 1:N
    if strcmpi(X.age_type(i),'14C') == 0       %non-14C age
        %make a truncated normal distribution
        norm_dist = makedist('Normal',X.age_common(i),X.age_err(i));
        if strcmpi(age_scale,'BC/AD') == 1 
            truncated_norm = truncate(norm_dist,A,B);
        else  
            truncated_norm = truncate(norm_dist,B,A);
        end    
        %calculate the pdf of the truncated normal distribution
        prob = pdf(truncated_norm,T);
    elseif strcmpi(X.age_type(i),'14C') == 1   %14C age 
        % retrieve the calibration dataset
        curve_cal_age = CalCurves.(char(X.cal_curve(i)))(:,1);
        curve_F14_val = CalCurves.(char(X.cal_curve(i)))(:,2);
        curve_F14_err = CalCurves.(char(X.cal_curve(i)))(:,3);
        hi_curve_cal_age = curve_cal_age(1):dt:curve_cal_age(end);
        hi_curve_cal_age = hi_curve_cal_age';
        hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
        hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
        % calculate the highest posterior probabilities
        a = (X.F14C_val(i) - hi_curve_F14_val).^2;
        b = 2*(X.F14C_err(i)^2 + hi_curve_F14_err.^2);
        c = sqrt(X.F14C_err(i)^2 + hi_curve_F14_err.^2);
        prob = exp(-a./b)./c;
        prob = flipud(prob);
    end
    % normalize the probabilities to 1
    prob = prob/sum(prob);
    cal_prob(:,i+1) = prob;
end
return;