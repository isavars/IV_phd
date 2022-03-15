function [resPkFr res1stPkFr resTM resAC]=spk_intrfreqautocorr(trialData,varargin)
% Instrinsic frequency analysis based on temporal autocorrelogram.
%
%       [intrFr intrFr_1stPk thetaMod AC]=spk_intrfreqautocorr(trialData)
%
% Outputs have fields .cells(1:nCell) (the values for each cell) and .eeg, (the values for the EEG).
%
% intrFr        - intrinsic frequency based on the FFT of the 0-500ms auto-correlogram.
% intrFr_1stPk  - intrinsic frequency based the frequency of the first peak of the 0-500ms auto-correlogram.
% thetaMod      - theta modulation, ratio of power of FFT of auto-correlogram in theta band, versus everything else.
% AC            - the autocorrelogram.
%
% Input options are specified with name - parameter comma-separated list.
% Parameters (+defaults) as follows:
%     opt.speedFilt=[0 400]; % In cm/s
%     opt.chunkLength=0.5;
%     opt.eegCh=1;  % EEG channel
%     opt.acBin = 5; % In ms
%     opt.acWin = 500; % In ms
%     opt.thetaBand = [4 10]; % Used to select the initial theta peak
%     opt.peakBand = 0.5; % Theta power is assesed +- peakBand Hz of the peak frequency
%     opt.smoothPS = 1.5; % In Hz - for smoothing power spectra
%     opt.smoothAC = []; % In ms
%     opt.lowCutOff=2; % Disregard frequencies below this when assesing modulation power.
%     opt.fig=0;    % Draw a figure with results
%     opt.hAxes=[]; % Pass handle to existing figure to draw into.
%     opt.cellIndex=1:length(trialData.cells);    % Can just look at particular cells to save time.
%     opt.analyseEEG = 1;                         % Can also not analyse EEG to save time.

opt.speedFilt=[0 400]; % In cm/s
opt.chunkLength=0.5;
opt.eegCh=1;  % EEG channel
opt.acBin = 2; % In ms
opt.acWin = 500; % In ms
opt.thetaBand = [6 12]; % Used to select the initial theta peak
opt.peakBand = 0.5; % Theta power is assesed +- peakBand Hz of the peak frequency
opt.smoothPS = 1.5; % In Hz - for smoothing power spectra
opt.smoothAC = []; % In ms
opt.lowCutOff=2; % Disregard frequencies below this when assesing modulation power.
opt.fig=1;
opt.hAxes=[];
opt.cellIndex=1:length(trialData.cells);    % Can just look at particular cells to save time. Should be numerical index, not logical index.
opt.analyseEEG = 1;                         % Can also not analyse EEG to save time.
for ii=1:2:length(varargin)
    opt.(varargin{ii}) = varargin{ii+1};
end
if opt.fig
    if isempty(opt.hAxes)
        hFig=gra_multiplot(length(opt.cellIndex),3,'plotsize',[6 4]); ud=get(hFig,'userdata'); ax=ud.axes_handle_array';
        set(hFig,'name',trialData.trialname);
    else
        ax=opt.hAxes;
    end
end

% Filter by speed (if requested) %
posSampRate=trialData.sample_rate;
[posInd,nPosSamp,chunkInd]=pos_speedfilter(trialData,'min',opt.speedFilt(1),'max',opt.speedFilt(2),'chunk',opt.chunkLength);
if nPosSamp<1 % Less than 1second of data
    resPkFr.cells=NaN(1,length(trialData.cells));   resTM.cells=resPkFr.cells;
    resPkFr.eeg=NaN;   resTM.eeg=NaN;
    return
end

resPkFr.cells=nan(size(opt.cellIndex));
res1stPkFr.cells=nan(size(opt.cellIndex));
resTM.cells=nan(size(opt.cellIndex));
resAC(1:length(opt.cellIndex))=deal({nan});

%%% Calculate chunked AC and Powerspec for each cell %%%
for ii=1:length(opt.cellIndex)
    cellData=trialData.cells(opt.cellIndex(ii));
    if isempty(cellData.st);    continue;    end;
    % Construct spike train histogram for autocorrelogrgams %
    spkTrHist=histc(cellData.st, 0:(opt.acBin/1000):trialData.dur );
    % Find the AC for each chunk %
    acChunk=NaN(length(unique(chunkInd)),opt.acWin/opt.acBin);
    for jj=1:length(unique(chunkInd))
        posIndChunk=posInd(chunkInd==jj);
        % Find Spike train for this chunk and find AC %
        %%% NB. Caswell had a good point here. Why don't we ignore those speed chunks with less than a 
        %%%     certain number of spikes?
        spkTrHistChunk=spkTrHist( ceil( (posIndChunk(1)/posSampRate)*(1000/opt.acBin) : (posIndChunk(end)/posSampRate)*(1000/opt.acBin) ) );
        [ac lags]=xcorr(spkTrHistChunk,'unbiased');
        winInd= lags>0 & lags<=(opt.acWin/opt.acBin);
        acChunk(jj,1:sum(winInd))=ac(winInd); % Select only the relevant short window
        acChunk(jj,:)=acChunk(jj,:).*length(posIndChunk);  % To weight by chunk length
    end
    % Mean AC over all chunks, then find Power Spectrum %
    ac=nansum(acChunk,1)./length(posInd);    % To weight by chunk length
    ac = ac - nanmean(ac);    % Mean normalise to reduce low-frequency components in power spectrum.
    [ps,freqs] = powerSpectrum(ac,2^16,(1000/opt.acBin),125);
    ps = smoothPowerSpec(ps,freqs,opt.smoothPS);
    ac=smoothPowerSpec(ac,opt.acBin:opt.acBin:opt.acWin,opt.smoothAC);
    resAC{ii}=ac;
    % Quantify theta modulation of Cell %
    [peakFreq thetaMod peakPower] = quantifyThetaProps(ps,freqs,opt.thetaBand,opt.peakBand,opt.lowCutOff);
    resPkFr.cells(ii)=peakFreq;  resTM.cells(ii)=thetaMod;
    % Get Estimate of frequency from first peak of AC - highest local maxima in theta band %
    lagTimes=opt.acBin:opt.acBin:opt.acWin;
    lagFreqs=(1./lagTimes).*1000; % In Hz
    [pk,pkInd]=findpeaks(ac);
    thetaPksLInd=( lagFreqs(pkInd)<opt.thetaBand(2) & lagFreqs(pkInd)>opt.thetaBand(1) );
    pkInd=pkInd(thetaPksLInd);   pk=pk(thetaPksLInd);
    [dum,maxPkInd]=max(pk);
    firstPkTime=lagTimes(pkInd(maxPkInd));
    firstPkFreq=lagFreqs(pkInd(maxPkInd));
    res1stPkFr.cells(ii)=firstPkFreq;
    % Plot Cell AC and PS %
    if opt.fig
        % AC figure %
        plot(ax(1,ii), opt.acBin:opt.acBin:opt.acWin,ac,'k-');  set(ax(1,ii),'ylim',[min(ac) max(ac)],'xlim',[0 opt.acWin]);
        hold(ax(1,ii),'on');  
        plot(ax(1,ii),[firstPkTime firstPkTime],get(ax(1,ii),'ylim'),'k-');
        text(firstPkTime,max(ac),{['1stPk=' num2str(firstPkFreq,'%3.1f')];['(' num2str(firstPkTime,'%d') 'ms)']},'verticalalignment','top','parent',ax(1,ii),...
            'fontunits','normalized','fontsize',0.1);
        hold(ax(1,ii),'off');  
        % PS Figure %
        plot(ax(2,ii), freqs(freqs<25),ps(freqs<25),'b-');      set(ax(2,ii),'xlim',[0 25]);
        if peakPower~=0; set(ax(2,ii),'ylim',[0 peakPower*1.5]); end
        hold(ax(2,ii),'on');   
        plot(ax(2,ii),[peakFreq peakFreq],get(gca,'ylim'),'k-');
        text(11,peakPower*1.5,{['PF=' num2str(peakFreq,'%3.1f')];['TM=' num2str(thetaMod,'%4.1f')]},'verticalalignment','top','parent',ax(2,ii),...
            'fontunits','normalized','fontsize',0.1);
        hold(ax(2,ii),'off');
    end
end
%%% Find the chunked powerspec for the EEG %%%
if ~isempty(trialData.eeg) && opt.analyseEEG
    eegChunk=repmat(NaN,1,length(trialData.eeg(opt.eegCh).eeg));
    eegChunkCount=1;
    for ii=1:length(unique(chunkInd))
        posIndChunk=posInd(chunkInd==ii);
        eegInd = repmat(posIndChunk.*5, 1, 5) - repmat(0:4, length(posIndChunk), 1); % Convert pos filter index
        eegInd = sort(reshape(eegInd, [], 1));                                       % .. to EEG filter index
        eegChunk( eegChunkCount:eegChunkCount+length(eegInd)-1 ) = trialData.eeg(opt.eegCh).eeg(eegInd);
        eegChunkCount=eegChunkCount+length(eegInd);
    end
    % EEG Powerspec %
    eegChunk=eegChunk(~isnan(eegChunk));
    [psEEG,freqsEEG]=powerSpectrum(eegChunk,2^19,posSampRate*5,25);
    psEEG = smoothPowerSpec(psEEG,freqsEEG,opt.smoothPS);
    [peakFreqEEG thetaModEEG] = quantifyThetaProps(psEEG,freqsEEG,opt.thetaBand,opt.peakBand,opt.lowCutOff);
    resPkFr.eeg=peakFreqEEG;   resTM.eeg=thetaModEEG; %  resPS.eeg=psEEG;
else
    resPkFr.eeg=NaN;   resTM.eeg=NaN;
end
% Plot EEG PS %
if opt.fig && opt.analyseEEG
    plot(ax(3,1),freqsEEG,psEEG,'k-');   set(ax(3,1),'xlim',[0 25]);
    hold(ax(3,1),'on');   
    plot(ax(3,1),[peakFreqEEG peakFreqEEG],get(ax(3,1),'ylim'),'k-');
    yLim=get(ax(3,1),'ylim');
    text(11,yLim(2),{['PF=' num2str(peakFreqEEG,'%3.1f')]; trialData.trialname},'verticalalignment','top','parent',ax(3,1),'fontunits','normalized','fontsize',0.15);
    hold(ax(3,1),'off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PS,freqs]=powerSpectrum(sigData,pad,sampRate,cutOffFr)
% Find power spectrum, of AC or EEG signal %
fTrans = fft(sigData, pad);
PS = fTrans .* conj(fTrans) / length(fTrans);
freqs = (0:(length(fTrans)-1)) * (sampRate/length(fTrans));
PS = PS(freqs<=cutOffFr);
freqs = freqs(freqs<=cutOffFr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [smPS]=smoothPowerSpec(power,freqs,siz)
% Smooth powerspec %
if isempty(siz);  smPS=power;  return;  end  % siz=[] means no smoothing to be done.
kernelLength = sum(freqs<=siz);
x = ((-siz/2)+(siz/kernelLength)) : (siz/kernelLength):(siz/2);
sigma = siz/4;
kernel = 10*exp(-0.5*(x/sigma).^2)';
smPS =  imfilter(power', kernel, 'symmetric')/sum(kernel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peakFreq,thetaMod,peakPower]=quantifyThetaProps(PS,freqs,thetaBand,peakBand,lowCutOff)
%%% Find peak frequency in theta band, and degree of modulation in this band %%%
% Low cut-off frequencies should be disregarded from all calculations, so get rid of them here.
PS=PS(freqs>lowCutOff);         
freqs=freqs(freqs>lowCutOff);
% To find theta peak, set all values outside theta band to 0 - easier to index peak frequency like this.
psTemp=PS;
psTemp( freqs<thetaBand(1)+1 | freqs>thetaBand(2)+1 ) = 0; % First, allow an extra 1Hz either side of window (as we still have to find local maxima)
localMaxLogInd=false(size(PS));                               % Find local maxima
[dum localMaxNumInd]=findpeaks(psTemp,'minpeakheight',0);     % Looking for >0 only saves time                          
localMaxLogInd(localMaxNumInd)=true;                          %
psTemp( freqs<thetaBand(1) | freqs>thetaBand(2) ) = 0;     % Now set ALL values outside window to 0.
% Find peak frequency and quantify %
psTemp(~localMaxLogInd)=0;     % Peak value must be a local maxima
[peakPower peakFreqInd]=max(psTemp);
peakFreq=freqs(peakFreqInd);
peakBandInd=( freqs>=(peakFreq-peakBand) & freqs<=(peakFreq+peakBand) );
warning('off', 'MATLAB:divideByZero');
thetaMod = mean(PS(peakBandInd)) / mean(PS(~peakBandInd));
warning('on', 'MATLAB:divideByZero');


