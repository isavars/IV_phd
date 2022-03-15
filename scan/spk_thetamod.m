function [thetaMax] = spk_thetamod(spkTr,trialLength)
% Estimate theta-modulation of cell from power spectrum of
% auto-correlogram.
%
%       [] = spk_thetamod(spkTr,thetaBand);

% data = evalin('base',SD.selData{1});
% spkTr = data.trials(SD.selTrial).cells(SD.selCell).st;
hFig = gra_multiplot(1,2,'plotsize',[6 4]);   ud=get(hFig,'userdata');

acWindow=0.5; % In seconds.
acDiv=250;
thetaBand = [4 10];
lowFreqCutOff = 2;

% New way of doing AC %
acBin=2;                                  % In ms.
spkTrHistInd=ceil((spkTr.*1000)/acBin);
spkTrHist=repmat(0,1,(trialLength.*1000)./acBin);
spkTrHist(spkTrHistInd)=1;
[ac lags]=xcorr(spkTrHist,'unbiased');
winInd=(lags>0 & lags<=(500/acBin));    % Use only 500ms window for FFT
acWin=ac(winInd);
lagsWin=lags(winInd);
axes(ud.axes_handle_array(1)); plot(lagsWin,acWin,'k-');

% spkAc = spk_crosscorr(spkTr,spkTr,acWindow*1000,'divisions',acDiv,'drawPlot',0,'hAxis',[]);%ud.axes_handle_array(1));
% figure; plot(1:acDiv+1,spkAc(251:501),'b-'); hold off;

% FFT of autocorrelogram %
spkFFT = fft(acWin, 1024);
pwrSpec = spkFFT .* conj(spkFFT) / length(spkFFT);
freqs = (0:(length(spkFFT)-1)) * ( (1000/acBin) / length(spkFFT) );
freqs = freqs(1:length(pwrSpec)/2);
pwrSpec = pwrSpec(1:length(pwrSpec)/2);

% Quantify %
thBandInd = freqs>thetaBand(1) & freqs<thetaBand(2);
nonThBandInd = freqs>lowFreqCutOff;
thetaRatio = mean(pwrSpec(thBandInd))  /  mean(pwrSpec( setdiff(find(nonThBandInd),find(thBandInd)) ));
thetaMax = max(pwrSpec(thBandInd));

% Plot %
axes(ud.axes_handle_array(2));
plot(freqs,pwrSpec,'b-');
ylim=get(gca,'ylim');
text(15,ylim(2),num2str(thetaRatio),'verticalalignment','top');
hold on;   plot([thetaBand(1) thetaBand(1)], [ylim(1) ylim(2)], 'k-');   
           plot([thetaBand(2) thetaBand(2)], [ylim(1) ylim(2)], 'k-');
           plot([lowFreqCutOff lowFreqCutOff], [ylim(1) ylim(2)], 'k-');    hold off;
set(gca,'xlim',[0 25]);










