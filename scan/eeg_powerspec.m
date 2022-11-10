function [peakFreq, peakPower, snr] = eeg_powerspec(eeg,varargin)
% Find power spectra of EEG data using FFT.
%
%   1) [peakFreq, peakPower, snr] = eeg_powerspec(eeg);       where eeg is a SCAn data .eeg struct.
%   2) [peakFreq, peakPower, snr] = eeg_powerspec(eegSamples, SR);   where eegSamples is eeg data, SR is sample rate.
%
% Power spectrum is smoothed using guassian window, default width 1.5Hz.
% Also finds and labels the highest local maxima in the theta band (default 4-10Hz).

if isstruct(eeg)
    eegSamp = ( double(eeg.eeg)./128 ) .* eeg.scalemax;
    sampRate = eeg.sample_rate;
elseif isnumeric(eeg)
    eegSamp = double(eeg);
    sampRate = varargin{1};
    varargin = varargin(2:end);
end
opt.meanNorm = 1;   % Mean normalise EEG before FFT.
opt.thetaBand = [5 11];
opt.snrBand = 1; % SNR of peak is calculated for band opt.snrBand width around peak.
opt.varBandWidth = 3;
opt.freqRes = [];   % Frequency resolution. If (eg) 0.01, final power spec will be at 0.01Hz resolution. If empty, resolution depends on signal length into fft.
opt.hfCutOff = 125;
opt.filterWidth = 2;
opt.hAxis = [];
opt.fig = 0;
opt.normToFreq = 0;
for ii=1:2:length(varargin)
    opt.(varargin{ii}) = varargin{ii+1};
end
% Allow for variable theta bands %
if length(opt.thetaBand)==1
    opt.thetaBand = [opt.thetaBand-(opt.varBandWidth/2) opt.thetaBand+(opt.varBandWidth/2)];
end

% Mean Normalise %
if opt.meanNorm
    eegSamp = eegSamp - mean(eegSamp,'omitnan');
end

% Remove NaNs that might be present %
eegSamp = eegSamp( ~isnan(eegSamp) );

% FFT %
% Second value should be a power of 2 greater in length than the eeg series, 2^19 is roughly the length of a 40 min trial.
fourierTrans = fft(eegSamp, 2.^(nextpow2(length(eegSamp)))); %result is complex
power = fourierTrans .* conj(fourierTrans) / length(fourierTrans);  % NB eqiv of taking the square modulus
% Frequencies at each power spectrum sample point %
freqs = ((0:(length(fourierTrans)-1)) * (sampRate/length(fourierTrans)))';
power = power(freqs<=opt.hfCutOff);
freqs = freqs(freqs<=opt.hfCutOff);

if opt.normToFreq
    power = power./freqs;
end

% Take mean power for each frequency at opt.freqRes sampling rate %
if ~isempty(opt.freqRes)
    dsChunkLength = sum( freqs<=opt.freqRes );
    dsNChunk      = floor( length(freqs) / dsChunkLength );
    dsInd         = bsxfun(@plus, (1:dsChunkLength)-1,  (1 : dsChunkLength : (dsChunkLength*dsNChunk-1))' );
    power         = nanmean( power(dsInd), 2 );
    freqs         = nanmean( freqs(dsInd), 2 );
end


% Smooth power %
if ~isempty(opt.filterWidth)
    kLength = sum( freqs<=opt.filterWidth );
    power   =  imfilter(power, fspecial( 'gaussian', [kLength 1], kLength/3 ), 'replicate');
end



% Find peak frequency in theta band, and SNR of peak %
thetaBandInd             = freqs>opt.thetaBand(1) & freqs<opt.thetaBand(2);
[peakPower, peakFreqInd] = max(power(thetaBandInd));
freqsTemp                = freqs(thetaBandInd);
peakFreq                 = freqsTemp(peakFreqInd);
% If peak is on the 'shoulder', find instead a local maximum within the thetaBand range (check first to save time with FINDPEAKS)
if peakFreqInd==1 || peakFreqInd==sum(thetaBandInd)
    [~,localMaxNumInd] = findpeaks(power,'minpeakheight',0);     % Looking for >0 only saves time
    localMaxLogInd     = false(size(power));  
    localMaxLogInd(localMaxNumInd) = true;
    localMaxThetaInd   = localMaxLogInd & thetaBandInd;
    if any( localMaxThetaInd )  % Need to check that there *is* actually a local max in the theta band - if not, stick with the 'shoulder' estimate we have.
        psTemp                     = power;  
        psTemp(~localMaxThetaInd)  = 0;
        [peakPower, peakFreqInd]   = max(psTemp);
        peakFreq                   = freqs(peakFreqInd);
    else
%         [peakPower, peakFreq, snr] = deal(nan); %% LM uncommented temporarily
    end
end
% Calculate SNR (only if there IS a local maxima in the theta band - see above for test) %
if ~isnan( peakFreq )
    peakBandInd = ( freqs>=(peakFreq-(opt.snrBand/2)) & freqs<=(peakFreq+(opt.snrBand/2)) );
    snr         = mean(power(peakBandInd)) / mean(power(~peakBandInd));
end


% % Plot power spectra %
% if opt.fig
%     if isempty(opt.hAxis);   figure;   else   axes(opt.hAxis);   end
%     plot(freqs,power,'b-');
%     fontspec = {'fontunits', 'normalized', 'fontsize', 0.075};
%     set(gca, 'ylim', [-peakPower*0.2 peakPower*1.3], 'xlim', [0 opt.hfCutOff], fontspec{1:end});
%     text('string',num2str(peakFreq, '%3.1f'),'position', [3 peakPower*1.3],'verticalalignment','top',fontspec{1:end});
%     text('string',num2str(snr, '%3.1f'),'position', [opt.hfCutOff peakPower*1.3],'verticalalignment','top','horizontalalignment','right',fontspec{1:end});
% end
   