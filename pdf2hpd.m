function CAL_AGE = pdf2hpd(pdfs,delta,X,age_scale) 
%% function for calculating the highest posterior density regions of the modeled calendar ages
%INPUT 
% pdfs: empirical pdf of the modeled calendar ages
% delta: nearest years the modeled calendar ages to be rounded
% X: structure containing information about ages of the SLIPs
% age_scale: scale of year (BC/AD, BP, B2K)
%OUTPUT
% CAL_AGE: structure containing modeled calendar ages including the 68.2% and
% 95.4% highest posterior density regions, the fraction of enclosed area, and the median probability age 
%%
N = length(X.age); 
CAL_AGE = struct('Sample_ID',[],'Unmodeled_age',[],'P68_2_Credible_intervals',[],'P95_4_Credible_intervals',[],'Median_prob_age',[],'Age_scale',[]);
cal_age = pdfs(:,1);
prob = pdfs(:,2:end);
for j = 1:N
    CAL_AGE(j).Sample_ID = X.sample_ID{j};
    CAL_AGE(j).Unmodeled_age = [num2str(X.age(j)),char(177),num2str(X.age_err(j))];
    [CAL_AGE(j).P68_2_Credible_intervals,CAL_AGE(j).P95_4_Credible_intervals,CAL_AGE(j).Median_prob_age] = hpds(cal_age,prob(:,j),delta);
    CAL_AGE(j).Age_scale = age_scale;
end
return;
%%
function [p68_2,p95_4,median_age] = hpds(cal_age,prob,delta)
%% function for estimating the one- and two-sigma highest posterior density regions of modeled calendar ages  
%INPUT
% cal_age: modeled calendar ages
% prob: probability of the modeled calendar ages   
% delta: nearest years the modeled calendar ages to be rounded
%OUTPUT
% p68_2: calendar ages at the 68.2% posterior probability
% p95_4: calendar ages at the 95.4% posterior probability
% median_age: calibrated age at median probability
%%
% Please see manuscript for more detailed information:
% Lougheed, B.C. & Obrochta, S.P. (2016). MatCal: Open Source Bayesian
% 14C Age Calibration in Matlab. Journal of Open Research Software. 4(1),
% p.e42. DOI: http://doi.org/10.5334/jors.130
%%
%% Calculating the highest posterior density regions
cal_age = cal_age(:);
prob = prob(:);
calprob = [cal_age prob];
hpd = calprob(:,1:2);
hpd = sortrows(hpd, 2);
hpd(:,3) = cumsum(hpd(:,2));
%one-sigma highest posterior density regions
hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)), :);
hpd68_2 = sortrows(hpd68_2,1);
ind1 = find(diff(hpd68_2(:,1)) > delta);
if isempty(ind1) == 1
	p68_2(1,1) = hpd68_2(end,1);
	p68_2(1,2) = hpd68_2(1,1);
	%p68_2(1,3) = sum(hpd68_2(1:end,2));
    p68_2(1,3) = 1; %fraction of enclosed area in the pdf
else
	z = 0;
	for i = 1:length(ind1)
		z = z + 1;
		indy1(z,1) = ind1(i);
		z = z + 1;
		indy1(z,1) = ind1(i)+1;
	end
	indy1 = [ 1 ; indy1; length(hpd68_2(:,1)) ];
	z=0;
	for i = 2:2:length(indy1)
		z = z+1;
		p68_2(z,1) = hpd68_2(indy1(i),1);
		p68_2(z,2) = hpd68_2(indy1(i-1),1);
		p68_2(z,3) = sum(hpd68_2(indy1(i-1):indy1(i),2));
	end
	p68_2 = flipud(p68_2);
    p68_2(:,3) = p68_2(:,3)/sum(p68_2(:,3)); %fraction of each enclosed area in the pdf
end
%two-sigma highest posterior density regions
hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)), :);
hpd95_4 = sortrows(hpd95_4,1);
ind2 = find(diff(hpd95_4(:,1)) > delta);
if isempty(ind2) == 1
	p95_4(1,1) = hpd95_4(end,1);
	p95_4(1,2) = hpd95_4(1,1);
	%p95_4(1,3) = sum(hpd95_4(1:end,2));
    p95_4(1,3) = 1; %fraction of enclosed area in the pdf
else
	z = 0;
	for i = 1:length(ind2)
		z = z + 1;
		indy2(z,1) = ind2(i);
		z = z + 1;
		indy2(z,1) = ind2(i)+1;
	end
	indy2 = [ 1 ; indy2; length(hpd95_4(:,1)) ];
	z=0;
	for i = 2:2:length(indy2)
		z = z+1;
		p95_4(z,1) = hpd95_4(indy2(i),1);
		p95_4(z,2) = hpd95_4(indy2(i-1),1);
		p95_4(z,3) = sum(hpd95_4(indy2(i-1):indy2(i),2));
	end
	p95_4 = flipud(p95_4);
    p95_4(:,3) = p95_4(:,3)/sum(p95_4(:,3)); %fraction of enclosed area in the pdf
end
%% calculate age at median probability
[~, median_ind] = min(abs(cumsum(calprob(:,2))-0.5));
median_age = round(calprob(median_ind(1),1));
return;