function [S,phi] = initials(T,Y,age_scale,G)
%% function for initializing model parameters
%INPUT
%T: initial values of the calendar age sequence
%Y: structure containing altitudinal info about the SLIPs
%age_scale: scale of year to be reported (e.g. BC/AD, BP, or B2K)
%G: temporal grid points
%OUTPUT
%S: initial values of sea level
%phi: initial values of sea-level rate
%% obtain the initial values of S
%get observational data
DT = unique(abs(diff(G)));
m = Y.rsl;
v = Y.rsl_err;
y = zeros(length(m),1);
true = length(unique(y)) == length(y);
while (~true)
    for i = 1:numel(m)
        norm_dist = makedist('Normal',m(i),v(i)); 
        truncated_norm = truncate(norm_dist,m(i)-v(i),m(i)+v(i));
        y(i) = random(truncated_norm,1,1);
    end 
    true = length(unique(y)) == length(y);
end
[t,~,idx] = unique(T,'stable');
val = accumarray(idx,y,[],@mean);
% add the origin to SLIPs
if strcmpi(age_scale,'BC/AD') == 1
    if t(end) ~= 2000
        t0 = [t; 2000];
        y0 = [val; 0];
    else
        t0 = t;   
        y0 = val;
    end    
else
    if t(end) ~= 0
        t0 = [t; 0];
        y0 = [val; 0];
    else
        t0 = t;
        y0 = val;
    end    
end
% Brownian bridge
%calculate mean
mu = interp1(t0,y0,G,'linear','extrap');
%calcualte variance
vr2 = zeros(length(G),1);                
if strcmpi(age_scale,'BC/AD') == 1
    for i = 1:numel(G)
        if G(i)<t0(1)
           vr2(i) = (t0(1)-G(i))^2/(t0(1)-G(1));
        end
    end
    for i = 1:numel(G)
        if G(i)>t0(end)
           vr2(i) = (G(i)-t0(end))^2/(G(end)-t0(end));
        end
    end
    for j = 1:numel(t)-1
        for i = 1:numel(G)
            if G(i)>t0(j) && G(i)<t0(j+1)
               vr2(i) = (t0(j+1)-G(i))*(G(i)-t0(j))/(t0(j+1)-t0(j));
            end 
        end
    end    
else
    for i = 1:numel(G)
        if G(i)>t0(1)
           vr2(i) = (G(i)-t0(1))^2/(G(1)-t0(1));
        end
    end
    for i = 1:numel(G)
        if G(i)<t0(end)
           vr2(i) = (t0(end)-G(i))^2/(t0(end)-G(end));
        end
    end
    for j = 1:numel(t)-1
        for i = 1:numel(G)
            if G(i)<t0(j) && G(i)>t0(j+1)
               vr2(i) = (t0(j)-G(i))*(G(i)-t0(j+1))/(t0(j)-t0(j+1));
            end 
        end
    end    
end
%draw samples
vr = sqrt(2*DT*vr2/(abs(G(1)-G(end)))/(length(t0)-1));
vr0 = min(vr(vr~=0));
S = zeros(length(G),1);
for i = 1:numel(G)
    if vr(i) == 0
        vr(i) = vr0;
    end    
    norm_dist = makedist('Normal',mu(i),vr(i)); 
    truncated_norm = truncate(norm_dist,mu(i)-vr(i),mu(i)+vr(i));
    S(i) = random(truncated_norm,1,1);
end    
%force the rsl to be 0 at "present" per the definition
if strcmpi(age_scale,'BC/AD') == 1
    if G(end) == 2000
        S(end) = 0;
    end    
else
    if G(end) == 0
        S(end) = 0;
    end    
end
%% obtain the initial value of phi
N = length(G)-1;
phi = zeros(N+1,1);
for i = 1:N+1
    if i == 1
       phi(i) = (S(i+1)-S(i))/DT;
    elseif i > 1 && i < N+1
       phi(i) = (S(i+1)-S(i-1))/(2*DT);
    else 
       phi(i) = (S(i)-S(i-1))/DT;
    end
end
end