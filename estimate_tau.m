function [ MiMat,OptTau] = estimate_tau(x,MaxDelay,Bins )
% Summary of this function goes here
%   get the optimal \tau to reconstruct phase space based on mutual info
%   x is the scalar time series samples, MaxDelay is the maximum delay
%   value and Bins is the no. of bins for histogram (= 16 is fine). Make 
% sure that x is an N X 1 vector

%detrend the data and get basic info
x = detrend(x); %detrend
TimePoints = size(x,1); %no of time points
N = size(x,2); % no of channels (at the moment works for single-channel data only - coloumn vector)

% set few parameters
TauMax = MaxDelay;
CorrRange = TimePoints - 2*TauMax;
BinEdge = ceil(CorrRange/(Bins));
SymbArray = zeros (2*TauMax + 1, N, CorrRange);
DataArray = x';

%compute entropies and auto mi

for t = 1:1:(2*TauMax + 1)
    Array = DataArray(:, t:CorrRange+t-1);
    EdgesSort = sort(Array,'ascend');
    Edges = EdgesSort(1:BinEdge:end);
    NoOfBins = length(Edges);
    for i = 1:1:NoOfBins
        Temp (i,:) = logical(Array(1,:) >= Edges(1,i));
    end
    SymbArray(t,:,:) = sum(Temp',2);
end


JointEntropy = 0.0 ;
for t = 1:1:CorrRange
    gFunc(t) = t*log(t);
end
Hist2D = int32(zeros(NoOfBins, NoOfBins));
MiMat = double(zeros(2*TauMax + 1, 1));



for t=1:1:(2*TauMax + 1)
    Tau = t - TauMax;
    for k=1:1:CorrRange
        Indexi = SymbArray(TauMax+1, 1, k);
        Indexj = SymbArray(t, 1, k);
        Hist2D(Indexi, Indexj) = Hist2D(Indexi, Indexj)+1;
    end
    
    Jointentropy = 0.0;
    for l= 1:1:NoOfBins
        for m = 1:1:NoOfBins
            if Hist2D(l,m) == 0
                Jointentropy = Jointentropy - gFunc(1);
            else
                Jointentropy = Jointentropy - gFunc(Hist2D(l,m));
            end
            Hist2D(l,m) = 0;
        end
    end
    
    Jointentropy = Jointentropy/CorrRange;
    Jointentropy = Jointentropy + log(CorrRange);
    
    MI = 0.0;
    MI = 2.0 * log( NoOfBins) - Jointentropy;
    MI = MI / log( NoOfBins);
    MiMat(t,1) = MI;
end


MiMatFromZero = MiMat(TauMax+1:end,:); % get only the mi values from lag 0 to max delay

%plot
semilogy(0:1:TauMax,MiMatFromZero,'b.-','LineWidth',2);
xlabel('lags')
ylabel('Mutual information');

% define optimal tau as first local minimum
DataInv = 1.001*max(MiMatFromZero) - MiMatFromZero;
[Minima,MinIdx] = findpeaks(DataInv);
OptTau = MinIdx(1)-1;
hold on
line([OptTau OptTau], get(gca, 'ylim'),'Color','r','LineWidth',2); 


end

