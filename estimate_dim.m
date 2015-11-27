function [FnnStat] = estimate_dim(x,OptTau,MinDim,MaxDim,r )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 

x = detrend(x);
Sigma = std(x,1);
FnnStat = zeros(MaxDim - MinDim + 1,1);

for m = MinDim:1:MaxDim
    Embedding_m = compute_psv(x,OptTau,m);
    Embedding_m1 = compute_psv(x,OptTau,m+1);
    Tree = KDTreeSearcher(Embedding_m(1:end-OptTau,:));
    [IDX,D] = knnsearch(Tree,Embedding_m(1:end-OptTau,:),'K',2,'distance','chebychev');
    D_m = D(:,2);
    IDX_m = IDX(:,2);
    
    D_m1 = max(abs(Embedding_m1 - Embedding_m1(IDX_m,:)),[],2);
    a = logical((D_m1 ./ D_m) > r);
    b = logical(D_m < (Sigma / r));
    
    FnnStat(m - MinDim+1,1) = sum(a.*b) / sum(b);
end

%plot
plot(MinDim:MaxDim,FnnStat,'b.-','LineWidth',2);
xlabel('dimension')
ylabel('FNN statistic');


end

