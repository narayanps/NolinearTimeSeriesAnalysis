function [C,T,L,R] = MovingWindowRNA(Data,ParamWin,Tau,Dim,Method)
%UNTITLED2 compute moving window network parameters for univariate data
%samples
% INPUT : Data - single data vector (1 X N)
%         ParamWin : structure comprising of
                     % WinLen in seconds - length of each window
                     % 
                     % OverLap - % overlap - number between 0 and 1
                     % SampFreq - Sampling frequency in Hetrz
%         Tau - Optimal lag  
%         Dim - Optimal dimension 
%         Method - string : ANN or FixedRR 
% OUTPUT : C,T,L,R - graph theoretic metrics - use BCT toolbox

ParamWin.WinStep = ParamWin.WinLen-(ParamWin.OverLap*ParamWin.WinLen);
[WindowedData,WinIndex] = MakeWindows(Data, ParamWin );
C = double(zeros(size(WindowedData,3),1));
T=C;
L=C;
R=C;

for i = 1:1:size(WindowedData,3)
    x = squeeze(WindowedData(1,:,i));
    x = detrend(x);
    x = double(x);
    [psa] = compute_psv(x',Tau,Dim);
    if (strcmp(Method,'ANN'))
           [A] = ComputeRecurrenceNetwork_ANN(psa,'sup_norm',10);
    elseif (strcmp(Method,'FixedRR'))
          [A,thr] = ComputeRecurrenceNetwork_fixedRR(psa,'sup_norm',0.02);
    end
    
    T(i,1) = transitivity_bu(A);
    C(i,1) = mean(clustering_coef_bu(A));
    R(i,1) = assortativity_bin(A,0);
    d = distance_bin(A);
    [L(i,1),efficiency,ecc,radius,diameter] = charpath(d);

end


end

