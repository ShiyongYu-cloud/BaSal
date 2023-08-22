function [] = plot_slips1(T_pdfs,CAL_AGE,Y,A,B,age_scale)
%% function for plotting SLIPs against modeled ages
%INPUT
% T_pdfs: array containing the empirical pdf of modeled calendar ages 
% CAL_AGE: highest posterior density regions of the modeled calendar ages
% Y: structure containing the altitudinal information of SLIPs
% A: early boundary of the study period
% B: late boundary of the study period
% age_scale: scale of the calibrated ages (BC/AD, BP, B2K)
%%
%% plot the SLIPs
cal_age = T_pdfs(:,1); 
prob = T_pdfs(:,2:end);
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
end