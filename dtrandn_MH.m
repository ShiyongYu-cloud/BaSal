function X = dtrandn_MH(X,Mu,Sigma,Mum,Mup)
%% function for sampling a truncated normal distribution bounded between Mum and Mup with 
% mean Mu and standard deviation Sigma 
% Modified after Dobigeon et al., IEEE TRANSACTIONS ON SIGNAL PROCESSING, 
% VOL. 57, NO. 11, NOVEMBER 2009
%--------------------------------------------------------------------------
%INPUT
% X: initial value of a variabe to be sampled
% Mu: mean of a variabe to be sampled
% Sigma: standard deviation of a variabe to be sampled
% Mum: lower bound of the variable
% Mup: upper bound of the variable
%OUTPUT
% X: updated value of the random variable
%%
Mu_new = Mu - Mum;
Mup_new = Mup -Mum;
if Mu<Mup
    Z = randnt(Mu_new,Sigma,1);
else
    delta = Mu_new - Mup_new;
    Mu_new = -delta;
    Z = randnt(Mu_new,Sigma,1);
    Z = -(Z-Mup_new );
end
    Z = Z+Mum;
    cond = (Z<=Mup) && (Z>=Mum);
    X = (Z.*cond + X.*(~cond));
return;
%%
function x = randnt(m,s,N)
% Generate N random numbers from a positive normal distribution with 
% mean M and standard deviation S
if s<0,  error('Standard deviation must be positive.'); end
if N<=0, error('N is wrong.'); end
Tindcand = [];
x = [];     
NN = N;
% Intersections
A  = 1.136717791056118;
mA = (1-A^2)/A*s;
mC = s * sqrt(pi/2);
while length(x)<NN
	if m < mA      % 4. Exponential distribution
	   a = (-m + sqrt(m^2+4*s^2)) / 2 / s^2;
       z = -log(1-rand(N,1))/a;
       rho = exp( -(z-m).^2/2/s^2 - a*(m-z+a*s^2/2) );
	elseif m <= 0  % 3. Normal distribution truncated at the mean
           z = abs(randn(N,1))*s + m;
           rho = (z>=0);
	elseif m < mC  % 2. Normal distribution coupled with the uniform one
           r = (rand(N,1) < m/(m+sqrt(pi/2)*s));
           u = rand(N,1)*m;
           g = abs(randn(N,1)*s) + m;
           z = r.*u + (1-r).*g;
           rho = r.*exp(-(z-m).^2/2/s^2) + (1-r).*ones(N,1);
    else           % 1. Normal distribution
        z = randn(N,1)*s + m;
        rho = (z>=0);
    end
	% Accept/reject the propositions
	reject = (rand(N,1) > rho);
	z(reject) = [];
	if ~isempty(z), x = [x ; z]; end
    N = N-length(z);
end
return;