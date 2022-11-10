function [stateInds] = getBrainStateHardThr( speed, eeg, varargin )
% Determine the brain state (run/sws/rem) on the basis of theta/delta ratio and speed.
% Unlike Zugaro's FMA toolbox, this function uses a set of hard thresholds to try and classify these
% (see below for details of how these are specified)
%
%       [stateInds] = getBrainStateHardThr( eeg, speed, prms )
%
% Output: stateInds.run, stateInds.sws, stateInds.rem are logical indexes, whose timebase is prms.stateStep.
%
% Analysis parameters:
%
% % Time windows for analysis %
% prms.posSR       = 50;          % Sample rates: I think these are completely standard, but there is the option ..
% prms.eegSR       = 250;         %    ..  to override these, if the need came about.
% prms.stateWindow = 1.6;         % States calculated in sliding windows this wide ..
% prms.stateStep   = 0.8;         % .. that slide by this much. 1.6 and 0.8 (sec) is the Csicsvari/Dupret default.
% prms.slidingSpeedWins = 0;      % Are speed windows sliding? (can make them not, possibly fair as speed already smoothed).
% % Parameters for spectrogram and frequency ratios %
% prms.thetaFr     = [];          % Theta frequency to use for theta/delta ratio - can be caller-supplied, but if empty, this function calculates using FFT.
% prms.thetaBW     = 1.5;         % Band width in which theta power is calculated: centred on prms.thetaFr.
% prms.deltaFr     = [];          %
% prms.deltaBW     = 2;
% prms.nTapers     = 2;           % 
% % Thresholds for state definitions %
% prms.speedThrRun   = 2.5;     % For run, must have speed faster than this ..
% prms.TDRatioThrRun = 0.5;     %  .. and TD ratio greater than this (makes sure we are in a 'theta' epoch).
% prms.speedThrSleep = 1;       % For both sleep states, speed must be less than this
% prms.TDRatioThrSWS = 0.75;    % SWS TD ratio must be less than this
% prms.TDRatioThrREM = 1;       % REM TD ratio must be greater than this.
% % Other %
% prms.brainStatePlots = 1;       % Make plots, for testing.
%
% %%% TODO - do we need a lower time limit for how long a rat can be in a state?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Parameters %
% Time windows for analysis %
prms.posSR       = 50;          % Sample rates: I think these are completely standard, but there is the option ..
prms.eegSR       = 250;         %    ..  to override these, if the need came about.
prms.stateWindow = 1.6;         % States calculated in sliding windows this wide ..
prms.stateStep   = 0.8;         % .. that slide by this much. 1.6 and 0.8 (sec) is the Csicsvari/Dupret default.
prms.slidingSpeedWins = 0;      % Are speed windows sliding? (can make them not, possibly fair as speed already smoothed).
prms.minEpochDur = 1.6;         % Brain state is only defined for states that are consistent for at least this long.
% Parameters for spectrogram and frequency ratios %
prms.thetaFr     = [];          % Theta frequency to use for theta/delta ratio - can be caller-supplied, but if empty, this function calculates using FFT.
prms.thetaBW     = 1.5;         % Band width in which theta power is calculated: centred on prms.thetaFr.
prms.deltaFr     = [];          %
prms.deltaBW     = 2;
prms.nTapers     = 2;           % 
% Thresholds for state definitions %
prms.speedThrRun   = 2.5;     % For run, must have speed faster than this ..
prms.TDRatioThrRun = 2;     %  .. and TD ratio greater than this (makes sure we are in a 'theta' epoch).
prms.speedThrREM   = 1;       % For sleep states, speed must be less than these ..
prms.speedThrSWS   = 2.5;     % .. REM needs a stricter threshold, otherwise RUN easily misclassified as REM.
prms.TDRatioThrSWS = 2;    % SWS TD ratio must be less than this
prms.TDRatioThrREM = 2;       % REM TD ratio must be greater than this.
% Other %
prms.brainStatePlots = 1;       % Make plots, for testing.
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Some preparation: get the number of windows:
trialDur = min( [length(eeg)./prms.eegSR, length(speed)./prms.posSR]  );  % Take the minimum of duration implied by speed and eeg, in case they don't match.
nWin     = floor( (trialDur-prms.stateWindow) / prms.stateStep );
% Get the theta frequency (if not supplied):
if isempty( prms.thetaFr )
    prms.thetaFr = eeg_powerspec(eeg, prms.eegSR, 'hfCutOff', 15, 'thetaBand',[5 11]);
end
% Get the delta frequency (if not supplied):
if isempty( prms.deltaFr )
    prms.deltaFr = eeg_powerspec(eeg, prms.eegSR, 'hfCutOff', 15, 'thetaBand',[1.5 4]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Get average speed per window %
if prms.slidingSpeedWins
    posWin = prms.stateWindow * prms.posSR; 
else
    posWin = prms.stateWindow/2 * prms.posSR;   % If stateWindow is /2, it means that we are not using overalapping (sliding) windows for speed.
end
posStep    = prms.stateStep * prms.posSR;
reshapeInd = bsxfun( @plus,  (posWin/2):posStep:(posStep*nWin),  (0:(posWin-1))' );
speedByWin = mean( speed(reshapeInd) ,  1  );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Get theta/delta per window %
eegWin           = prms.stateWindow * prms.eegSR;
eegStep          = prms.stateStep * prms.eegSR;
reshapeInd       = bsxfun( @plus,  (eegWin/2):eegStep:(eegStep*nWin),  (0:(eegWin-1))' );
[eegSpect, spFr] = pmtm(  eeg(reshapeInd),  prms.nTapers, 0.1:0.1:10, prms.eegSR);  % NOTE: frequencies are hard-coded here! (3rd argument).
thetaInd         = spFr>=(prms.thetaFr-prms.thetaBW/2) & spFr<=(prms.thetaFr+prms.thetaBW/2);
deltaInd         = spFr>=(prms.deltaFr-prms.deltaBW/2) & spFr<=(prms.deltaFr+prms.deltaBW/2);
TDRatio          = mean( eegSpect( thetaInd, : ), 1 ) ./ mean( eegSpect( deltaInd, : ), 1 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Set brain state output %
stateInds.run = speedByWin>=prms.speedThrRun & TDRatio>=prms.TDRatioThrRun;
stateInds.sws = speedByWin<prms.speedThrSWS & TDRatio<prms.TDRatioThrSWS;
stateInds.rem = speedByWin<prms.speedThrREM & TDRatio>prms.TDRatioThrREM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5) Filter state epochs to only use those of duration > prms.minEpochDur
f = {'run','sws','rem'};
if ~isempty( prms.minEpochDur )
    for ii=1:3
        stateDiff  = diff([0 double( stateInds.(f{ii}) ) 0]);
        epochStart = find(stateDiff == 1);
        epochEnd   = find(stateDiff == -1) - 1;
        epochDur   = (epochEnd - epochStart + 1) .* prms.stateStep;
        shortEpochList = find( epochDur < prms.minEpochDur );
        for jj=1:length(shortEpochList)
            stateInds.(f{ii})( epochStart(shortEpochList(jj)) : epochEnd(shortEpochList(jj)) ) = false;
        end
    end
end
% Add the TDRatio to output%
stateInds.TDRatio = TDRatio;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing plots %
if prms.brainStatePlots
    figure; 
    % IMPORTANT! each time bin is normalised here for plotting (otherwise difficult to see trends). %
    normSpect = bsxfun( @rdivide, eegSpect, sum( eegSpect,1 ) );
    
    % 1) Plot spectrogram %
    subplot(4,1,1); imagesc( normSpect );
    hold on;
    lineY = [prms.thetaFr-prms.thetaBW/2, prms.thetaFr+prms.thetaBW/2, prms.deltaFr-prms.deltaBW/2, prms.deltaFr+prms.deltaBW/2,];
    for ii=1:length(lineY);   plot( get(gca,'xlim'), [1 1].*lineY(ii).*10, 'w:' );   end

    % 2) Plot speed %
    subplot(4,1,2); plot( speedByWin );
    set( gca, 'xlim', [0 length(speedByWin)] );
    hold on;
    lineY = prms.speedThrSWS;
    for ii=1:length(lineY);   plot( get(gca,'xlim'), [1 1].*lineY(ii), 'k:' );   end
    
    % 3) Plot ratio %
    subplot(4,1,3); plot( TDRatio );
    set( gca, 'xlim', [0 length(speedByWin)], 'ylim', [0 10] );
    hold on;
    lineY = prms.TDRatioThrSWS;
    for ii=1:length(lineY);   plot( get(gca,'xlim'), [1 1].*lineY(ii), 'k:' );   end
    
    % 4) Plot state %
    stateVect = nan( size(stateInds.run) );
    stateVect( stateInds.run ) = 1;
    stateVect( stateInds.sws ) = 2;
    stateVect( stateInds.rem ) = 3;
    
    subplot(4,1,4); lh=plot( stateVect );
    set(gca, 'ylim', [0 4],'xlim', [0 length(speedByWin)] );
    set(lh,'color','k','linewidth',2);

%     figure;   plot(  speedByWin, TDRatio, 'bo'  );   set( gca, 'xlim', [0 15], 'ylim', [0 5] );
end

    