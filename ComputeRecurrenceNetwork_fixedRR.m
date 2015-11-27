function [A,thr] = ComputeRecurrenceNetwork_fixedRR(psv,dist_metric,rr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = size(psv,1);
distance = double(zeros(N,N));


%% compute the distance matrix%%
if strcmp(dist_metric,'euc_norm')
    
    for j = 1:1:N
        psv2 = repmat(psv(j,:),N,1);
        distance(:,j) = sqrt(sum((psv-psv2).^2,2));
    end
    
elseif strcmp(dist_metric,'man_norm')
    for j = 1:1:N
        psv2 = repmat(psv(j,:),N,1);
        distance(:,j) = sum(abs(psv-psv2),2);
    end
    
elseif strcmp(dist_metric,'sup_norm')
    for j = 1:1:N
        psv2 = repmat(psv(j,:),N,1);
        distance(j,:) = max(abs(psv-psv2)');
    end
end


%% construct the network
R = zeros(N,N);
tmp = reshape(distance,1, N*N);
tmp = sort(tmp);
thr = tmp(ceil(rr * length(tmp)));
R(distance < thr) = 1;
A = R;
A(logical(eye(size(A)))) = 0;
end






