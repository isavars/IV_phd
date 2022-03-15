function [H] = gra_plotwf(amps, scalemax, mode, varargin)
% Plot the average waveform of a cell. 
%
%       gra_plotwf(ampMatrix, scaleMax, 'max');
%       gra_plotwf(ampMatrix, scaleMax, 'all');
%
% .. to plot the maximum amplitude channel, or all channels, respectively.
%
% 'ampMatrix' can be either the full matrix of tetrode amp data (format=[nSpike, sample_1:50, tetChannel_1:4]),
% or alternatively can be already a matrix of meaned waveforms (format=[sample_1:50, tetChannel_1:4]).
%
% 'scalemax' is a 1x4 vector giving the max amplitude of the oscilloscope on each tetrode channel.
%
%       gra_plotwf(ampMatrix, scaleMax, 'mode', hAxis);
%
%  .. to pass a specific axis handle. For mode = 'all', hAxis should
% be a vector of 4 axes handles.

% Format for voltages is [spike(1:n), sample(1:50), channel(1:4)].
% Format for gains is {[gain_ch_1], [gain_ch_2], [gain_ch_3], [gain_ch_4]}.


%Note from LM: bug when trying to plot cells with 0 spikes
charWF = 'char';

if isempty(amps)
    % First, protect against no spike inputs %
    mean_wf = nan(50,4);
    maxCh = 1;
    charWF = 'dontChar';
else
    % otherwise, calculate mean and get max channel %
    if ndims(amps)==3  
        mean_wf = squeeze(mean(double(amps)));              % For matrices with every amp sample, take mean
    elseif ndims(amps)==2 && all(size(amps)==[50 4])
        mean_wf = double(amps);                             % In the case that amps is already the meaned waveforms, from scan.
    else
        error('Incorrect format for spike amplitude data.');
    end
    [dummy, maxCh] = max( max(mean_wf) - min(mean_wf) );
end
% Convert int8 values to actual voltages %
for ii=1:4
    mean_wf(:,ii) = (mean_wf(:,ii)./256) .* 2 .*scalemax(ii);
end
% Draw waveforms into axes. If axes are not supplied will be created here (either 1x1 or 1x4, depending on all or max channels) %
switch mode
    case 'max'
        if isempty(varargin)
            hAx=axes;
        else
            hAx=varargin{1};
        end
        H.handles = draw_wf_in_axes(mean_wf(:,maxCh), scalemax(maxCh), hAx, charWF);     % Last arg says whether or not to characterise WF
    case 'all'
        for ii=1:4
            if isempty(varargin)
                hAx=subplot(4,1,ii);
            else
                hAx=varargin{1}(ii);
            end
            if ii==maxCh
                H.handles(ii) = draw_wf_in_axes(mean_wf(:,ii), scalemax(ii), hAx,  charWF);
            else
                H.handles(ii) = draw_wf_in_axes(mean_wf(:,ii), scalemax(ii), hAx, 'dontChar' );
            end
        end
end

% Assign Output %
H.maxChannel = maxCh;

%%% This sub-function actually does the plotting %%%
function [handles] = draw_wf_in_axes(mean_wf, scalemax, hAx, characterise)
% Basic Plot %
plot(hAx, 1:50, mean_wf);
hold(hAx,'on');
plot(hAx,[1 50], [0 0], 'k-');
scalemax = ceil(scalemax);
set(hAx, 'xtick', [], 'xlim', [1 50], 'ylim', [-scalemax, scalemax], 'ytick',[]);
set(hAx, 'fontunits', 'normalized', 'fontsize', 0.1);
handles.yAxis_text = text(0.5,scalemax,[num2str(scalemax) 'uV'],'parent',hAx,'verticalalignment','bottom','horizontalalignment','right','fontunits','normalized','fontsize',0.125,'rotation',90);

% Characterise WF (only for max amp channel %
if strcmp(characterise,'char')
    % Find WF peak (highest amplitude local maxima above 0uV) %
    [peakAmpTemp,peakIndTemp] = findpeaks(mean_wf,'minpeakheight',0);
    [peakAmp,tempInd] = max(peakAmpTemp);
    peakInd = peakIndTemp(tempInd);
    % Find WF troughs (lowest local minima, we need to get two, both before and after the peak) %
    [preTrfAmpTemp,preTrfIndTemp] = findpeaks(-mean_wf(1:peakInd));
    if ~isempty(preTrfIndTemp)
        [preTrfAmp,tempInd] = min(-preTrfAmpTemp);
        preTrfInd = preTrfIndTemp(tempInd);
    else
        preTrfInd = NaN;   preTrfAmp = NaN;
    end
    [postTrfAmpTemp,postTrfIndTemp] = findpeaks(-mean_wf(peakInd:end));
    if ~isempty(postTrfIndTemp)
        [postTrfAmp,tempInd] = min(-postTrfAmpTemp);
        postTrfInd = postTrfIndTemp(tempInd) + peakInd - 1;
    else
        postTrfInd = NaN;   postTrfAmp = NaN;
    end
    
    % Print spike amplitude and spike width (peak-to post trough) on plot %
    str{1} = [num2str((max(mean_wf)-min(mean_wf)),'%3.0f') 'uV'];
    if ~isnan(postTrfInd)
        str{2} = [num2str(   (postTrfInd-peakInd)*0.0208*1000    , '%3.0f') 'us'];
    else
        str{2} = 'n/a';
    end
    handles.props_text = text('position',[1 0],'units','normalized','HorizontalAlignment','right','string',str,'FontUnits','normalized','VerticalAlignment','bottom','fontsize',0.1,'parent',hAx);
    
    % Plot peaks and troughs %
    handles.peak_line = plot(hAx,[peakInd peakInd],[-scalemax scalemax],'r:');
    handles.trough_line(1) = plot(hAx,[preTrfInd preTrfInd],[-scalemax scalemax],'m:');
    handles.trough_line(2) = plot(hAx,[postTrfInd postTrfInd],[-scalemax scalemax],'m:');
else
    handles.props_text = nan;
    handles.peak_line = nan;
    handles.trough_line(1) = nan;
    handles.trough_line(2) = nan;
end











        