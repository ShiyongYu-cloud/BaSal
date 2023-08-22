function [] = save_mcmc(T_samples,S_samples,phi_samples,X,G,age_scale)
%% function for saving the results from MCMC run to files
%INPUT
% T_samples: caldendar ages generated from MCMC run
% S_samples: RSL generated from MCMC run 
% phi_samples: rate of RSL changes generated from MCMC run
% X: structure containing information about the ages of the SLIPs
% G: temporal grid points
% age_scale: scale of year to be reported (BC/AD, BP, B2K)
%OUTPUT
% mcmc_T.txt: data file of calendar ages generated from MCMC run
% mcmc_S.txt: data file of RSL generated from MCMC run 
% mcmc_phi.rxt: data file of rate of RSL changes generated from MCMC run
%%
%% save T 
%creat a header
[M,N] = size(T_samples); 
header1 = cell(1,N+1); 
header1{1,1} = 'MCMC no.';
for j = 1:N
    header1{1,1+j} = strcat(X.sample_ID{j},'(',age_scale,')');
end
% write header and data to the file
fmt = repmat('%25d ',1,N+1); 
fmt = [fmt,'\n'];
fid1 = fopen('mcmc_T.txt','wt');
fprintf(fid1,'%25s ',header1{1:end});
fprintf(fid1,'\n');
for k = 1:M
    fprintf(fid1,fmt,[k T_samples(k,:)]);
end    
fclose(fid1);
%% save S
P = length(G);
%creat a header
header3 = cell(1,P+1); 
header3{1,1} = 'MCMC no.';
for j = 1:P
    header3{1,1+j} = strcat(num2str(G(j)),'(',age_scale,')');
end
% write header and data to the file
fmt = repmat('%25s ',1,P-1); 
fmt = ['%25d ' fmt '%25d'];
fid2 = fopen('mcmc_S.txt','wt');
fprintf(fid2,'%25s ',header3{1:end});
fprintf(fid2,'\n');
for k = 1:M
    fprintf(fid2,fmt,[k S_samples(k,:)]); 
    fprintf(fid2,'\n');
end    
fclose(fid2);
%% save phi
%creat a header
header3 = cell(1,P+1); 
header3{1,1} = 'MCMC no.';
for j = 1:P
    header3{1,1+j} = strcat(num2str(G(j)),'(',age_scale,')');
end
% write header and data to the file
fmt = repmat('%25s ',1,P-1); 
fmt = ['%25d ' fmt '%25d'];
fid3 = fopen('mcmc_phi.txt','wt');
fprintf(fid3,'%25s ',header3{1:end});
fprintf(fid3,'\n');
for k = 1:M
    fprintf(fid3,fmt,[k phi_samples(k,:)]); % padd 0 to the last grid point
    fprintf(fid3,'\n');
end    
fclose(fid3);
return;