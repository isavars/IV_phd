function [rtn] = spk_1stMomentAC(stA,acWin,acBin,trDur)
% Calculates first moment (i.e. mean value) of Auto-correlogram (in ms).
%
%       [rtn] = spk_1stMomentAC(spikeTrain,ACWin,ACBin,trialLength)
%
% spikeTrain should contain the list of spike times (in seconds), 
% ACBin and ACWin should be specified in ms, and trialLength in seconds.
%
% See Csicsvari J, Hirase H, Czurkó A, Mamiya A, Buzsáki G.1999 for reference values in HPC.

if length(stA) == 1
    rtn = NaN;
    return;
end    

% Make spike train vector %
spkTrHistA=histc(stA,0:(acBin/1000):trDur);
% AC %
[ac lags]=xcorr(spkTrHistA,spkTrHistA,(acWin/acBin),'unbiased');
% Select only the relevant short window %
ac(lags==0) = NaN;
winInd= lags>0 & lags<=(acWin/acBin);
ac=ac(winInd);
lags=lags(winInd);

%get first moment (i.e. mean in ms) of ac. 
%mean = sum(x(i-end)*p(x(i-end))). add binsize/2 to center over bins
rtn = ( lags+(acBin/2) ) * (ac / sum(ac));