function [WindowedData,WinIndex] = MakeWindows( EEG_Data, ParamWin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


WinLen = ParamWin.WinLen;
WinStep = ParamWin.WinStep;
SampFreq = ParamWin.SampFreq;
Srate = 1/SampFreq; %sampling rate

WinLen = WinLen - rem(WinLen,Srate); %compute window length to nearest multiple of srate
WinStep = WinStep - rem(WinStep,Srate); %compute window step to nearest multiple of srate
WinLenPnts  = round(WinLen*SampFreq); % window size in points
WinStepPnts = round(WinStep*SampFreq);  % windowstep in points

[Chans, Pnts] = size(EEG_Data);
WinStartIdx  = 1:WinStepPnts:(Pnts-WinLenPnts)+1;
NumWins   = length(WinStartIdx);
WinIndex = zeros(NumWins,2);


WindowedData = zeros(Chans,WinLenPnts,NumWins);

for t=1:NumWins
    Winpnts = WinStartIdx(t):WinStartIdx(t)+WinLenPnts-1;
    WinIndex(t,1) = Winpnts(1);
    WinIndex(t,2) = Winpnts(end);
    WindowedData(:,:,t) = EEG_Data(:,Winpnts);
end

end

