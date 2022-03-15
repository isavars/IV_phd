function [P] = spk_characterisewaveform(amps, scalemax)
% Characterise the waveform of a cell. 
%
%       PROPS = spk_characterisewaveform(ampMatrix, scaleMax);
%
% 'ampMatrix' can be either the full matrix of tetrode amp data (format=[nSpike, sample_1:50, tetChannel_1:4]),
% or alternatively can be already a matrix of meaned waveforms (format=[sample_1:50, tetChannel_1:4]).
%
% 'scalemax' is a 1x4 vector giving the max amplitude of the oscilloscope on each tetrode channel.
%
% P.max2MinAmp          % WF max to WF min amplitude, on maximum channel
% P.pk2TrAmp            % Peak to (post-peak) trough amplitude, on maximum channel
% P.pk2TrTime           % Peak to (post-peak) trough time (i.e. spike width)
% P.hhWidth             % Half-height width. Half-height refers to half of the amplitude of the peak above zero. 
% P.qhWidth             % Quarter height width (see half-height width)
% P.zcWidth             % Spike width between up- and down-going zero crossings.
% P.pre2PkAmp           % Pre-peak trough to peak amplitude
% P.pre2PkTime          % Pre-peak trough to peak time
% P.pre2PostAmpDiff     % Difference between pre-peak and post-peak trough
% P.xChanAmpRatio       % Ratio between amplitude on maximum channel, and next largest channel. Amplitudes on other channels (not max amp channel)
                        %    are calculated using the voltages at peak and trough times defined using the max amplitude channel. The amplitude can be
                        %    either pre-peak trough to peak, or peak to post-peak trough, depending on which is larger.


% Pre-assign PROPS output %
P.max2MinAmp = nan;          % WF max to WF min amplitude, on maximum channel
P.pk2TrAmp = nan;            % Peak to (post-peak) trough amplitude, on maximum channel
P.pk2TrTime = nan;           % Peak to (post-peak) trough time (i.e. spike width)
P.hhWidth = nan;             % Half-height width. Half-height refers to half of the amplitude of the peak above zero. 
P.qhWidth = nan;             % Quarter height width (see half-height width)
P.zcWidth = nan;             % Spike width between up- and down-going zero crossings.
P.pre2PkAmp = nan;           % Pre-peak trough to peak amplitude
P.pre2PkTime = nan;          % Pre-peak trough to peak time
P.pre2PostAmpDiff = nan;     % Difference between pre-peak and post-peak trough
P.xChanAmpRatio = nan;       % Ratio between amplitude on maximum channel, and next largest channel. Amplitudes on other channels (not max amp channel)
                             % are calculated using the voltages at peak and trough times defined using the max amplitude channel. The amplitude can be
                             % either pre-peak trough to peak, or peak to post-peak trough, depending on which is larger.




% Protect against no spike inputs %
if isempty(amps);   return;   end

% Calculate mean waveform %
if ndims(amps)==3 && size(amps,2)==50 && size(amps,3)==4
    mean_wf = squeeze(mean(double(amps)));              % For matrices with every amp sample, take mean
elseif ndims(amps)==2 && all(size(amps)==[50 4])
    mean_wf = double(amps);                             % In the case that amps is already the meaned waveforms, from scan.
else
    error('Incorrect format for spike amplitude data.');
end
for ii=1:4;   mean_wf(:,ii) = (mean_wf(:,ii)./256) .* 2 .*scalemax(ii);   end    % Convert int8 sample numbers to voltages.

% Get max2min amplitude, and the maximum amplitude channel %
[P.max2MinAmp, maxCh] = max( max(mean_wf) - min(mean_wf) );
maxAmpWF = mean_wf(:,maxCh);

%%%% Find waveform peaks and troughs on the max amplitude (max2minAmp) channel %%%%
% Find WF peak (highest amplitude local maxima above 0uV) %
[peakAmpTemp,peakIndTemp] = findpeaks(maxAmpWF,'minpeakheight',0);
[peakAmp,tempInd] = max(peakAmpTemp);
peakInd = peakIndTemp(tempInd);
% For very bad waveforms, can be no local max on largest channel. In this case, give up and return all as NaN. %
if isempty(peakAmp);    return;    end
% Find WF troughs (lowest local minima, we need to get two, both before and after the peak) %
[preTrfAmpTemp,preTrfIndTemp] = findpeaks(-maxAmpWF(1:peakInd));
if ~isempty(preTrfIndTemp)
    [preTrfAmp,tempInd] = min(-preTrfAmpTemp);
    preTrfInd = preTrfIndTemp(tempInd);
else
    preTrfInd = NaN;   preTrfAmp = NaN;    % Troughs can also not exist (if there is no local minima in the 1ms spike sample window %
end
[postTrfAmpTemp,postTrfIndTemp] = findpeaks(-maxAmpWF(peakInd:end));
if ~isempty(postTrfIndTemp)
    [postTrfAmp,tempInd] = min(-postTrfAmpTemp);
    postTrfInd = postTrfIndTemp(tempInd) + peakInd - 1;
else
    postTrfInd = NaN;   postTrfAmp = NaN;  % Troughs can also not exist (if there is no local minima in the 1ms spike sample window %
end

%%%%% Using these peaks and troughs, characterise the waveform %%%%%
% Peak-to-(post)trough amplitude and spike width %
P.pk2TrAmp = peakAmp - postTrfAmp;
P.pk2TrTime = (postTrfInd - peakInd)*0.0208;

% Pre-trough to peak amplitude and width %
P.pre2PkAmp = peakAmp - preTrfAmp;
P.pre2PkTime = (peakInd - preTrfInd)*0.0208;

% Pre-trough to post-trough amplitude %
P.pre2PostAmpDiff = preTrfAmp - postTrfAmp;

% Half-height width. Half-height refers to amplitude of peak above zero. Look for            %
% the crossings of half height between the pre-peak trough (or time=0, if not existing)      %
% and the peak, and between the peak and the post-peak trough (or time=1ms, if not existing) %
% ADDITION, TW, 07/11/12. Adapt this code  to also return 25% height spike width, and zero-crossing spike width. %
hh = peakAmp / 2;
qh = peakAmp / 4;   % Quarter height
% Up-going phase crossing %
if ~isnan(preTrfInd);   upPhaseStart = preTrfInd;   else   upPhaseStart=1;   end
upPhaseInd=false(1,50);   upPhaseInd(upPhaseStart:peakInd) = true;
upPhase = maxAmpWF;       upPhase(~upPhaseInd) = nan;
[dum,upCrossInd] = nanmin(abs(upPhase - hh));
[dum,upCrossIndQH] = nanmin(abs(upPhase - qh));
if min(upPhase) > 0;   upCrossIndZC = nan;   else   [dum,upCrossIndZC] = nanmin(abs(upPhase));   end  % Zero crossing width not defined if WF always +ve
% Down-going phase crossing %
if ~isnan(postTrfInd);   downPhaseEnd = postTrfInd;   else   downPhaseEnd=50;   end
downPhaseInd=false(1,50);   downPhaseInd(peakInd:downPhaseEnd) = true;
downPhase = maxAmpWF;       downPhase(~downPhaseInd) = nan;
[dum,downCrossInd] = nanmin(abs(downPhase - hh));
[dum,downCrossIndQH] = nanmin(abs(downPhase - qh));
if min(downPhase) > 0;   downCrossIndZC = nan;   else   [dum,downCrossIndZC] = nanmin(abs(downPhase));   end  % Zero crossing width not defined if WF always +ve
% Final half-height (and quarter height) %
P.hhWidth = (downCrossInd - upCrossInd) * 0.0208;
P.qhWidth = (downCrossIndQH - upCrossIndQH) * 0.0208;
P.zcWidth = (downCrossIndZC - upCrossIndZC) * 0.0208;


% Cross-channel amplitude ratio %
[amp,tempInd] = nanmax([P.pk2TrAmp P.pre2PkAmp]);  % Get largest of pre- and post- trough amplitudes
if ~isnan(amp)  % If there is a defined amplitude
    if tempInd==1                                  
        trfInd = postTrfInd;
        maxChAmp = P.pk2TrAmp;
    else
        trfInd = preTrfInd;
        maxChAmp = P.pre2PkAmp;
    end     
    channelAmps = mean_wf(peakInd,:) - mean_wf(trfInd,:);
    channelAmps(maxCh) = nan;
    nextBigChAmp = nanmax(channelAmps);
    P.xChanAmpRatio = nextBigChAmp ./ maxChAmp;
end









