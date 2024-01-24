%round A and B to the nearest 1000 years
function [] = plot_age_prior(A,B,X,Y,age_scale)
A = round(A/1000)*1000; 
B = round(B/1000)*1000;
%divid [A B] into N bins
N = length(X.age);
delta_t = round((A-B)/N);
cal_age = zeros(N+1,1);
prob = zeros(N+1,N);
for i = 1:N+1
    if strcmpi(age_scale,'BC/AD') == 1 
        cal_age(i,1) = A+(i-1)*delta_t;
    else
        cal_age(i,1) = A-(i-1)*delta_t;
    end
end
for i = 1:N
    prob(i:i+1,i) = 1/abs(A-B);
end
prob_u = zeros(N+1,N);
for i = 1:N
    factor = Y.rsl_err(i)/max(prob(:,i));
    prob_u(:,i) = Y.rsl(i) + prob(:,i)*factor; % blow up to the upper limit of error
    ind1 = prob_u(:,i) > Y.rsl(i);
    early = max(cal_age(ind1));
    late = min(cal_age(ind1));
    ind2 = (cal_age <= early) & (cal_age >= late);
    calage = cal_age(ind2);
    CALAGE = [calage(1); calage; calage(end); calage(1)];
    PROB_u = [Y.rsl(i); prob_u(ind2,i); Y.rsl(i); Y.rsl(i)];
    %plot(CALAGE,PROB_u);
    fill(CALAGE,PROB_u,[1 1 1]*0.85,'EdgeColor','none');  % plot pdfsend
    hold on
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
end