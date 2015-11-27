function [R] = ComputeRecurrenceNetwork_ANN(psv,dist_metric,K)
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

[A B] = sort(distance,2,'ascend');
sorted_neighbors = B(:,1:end);
N = size (distance,1);
R = zeros(N,N);
order = 1:1:N;

for i = 1:1:K
    for j = 1:1:N        
        l = order(j);
        k = i + 1;
        while (R(l,sorted_neighbors(l,k)) == 1 && k < N)
            k = k+1;
        end
        R(l,sorted_neighbors(l,k)) = 1;
        R(sorted_neighbors(l,k),l ) = 1;
    end
end


end






