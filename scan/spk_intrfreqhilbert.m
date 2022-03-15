function [spikePhase,intrinsicFreq] = spk_intrfreqhilbert(spikeTimes, instPhase250hz, cycle, validCycle, maxElapsedCycles, sampleRate)
%INTRINSIC_FREQ_HILBERT Calcs spike phase and intrinsic freq with respsect to theta LFP
%       [spikePhase,intrinsicFreq] = spk_intrfreqhilbert(spikeTimes, instPhase250hz, cycle, validCycle, maxElapsedCycles, sampleRate)
% Uses Hilbert transform of theta LFP to determine the phase of all spikes (theta cycles
% with low power [defined by validCycle] are set to nan). Then by comparing semi-adjacent
% theta cycles determine what intrinsic frequency the spikes are being emited at to get
% the observed change in phase. Comparison is only made for theta cyles that contain
% spikes are are equal to or less than maxElapsedCycles apart (three is probably a good
% value)

% ARGS
% spikeTimes - time of each spike in seconds (as reported by DACQ)
% instPhase250hz - [lengthEEG x 1] phase for each point in the eeg sampled at 250hz should
%               not be 'wrapped' i.e. phase should not be bounded between 0 and 2pi
% cycle - [lengthEEG x 1] vector of integeres indicating for each pt. in the EEG which
%           cycle it belongs to e.g [0,0,0,1,1,1,1,2 ...] cycles are divided at phase 0 to
%           2pi where 0 is the peak of the signal
% validCyle - list of cyles that full fill acceptance criteria (normally based on power)
%           e.g. [1,2,4,5,8, ..]
% maxElapsedCycles [x] - phase advance & intrinsic freq is computed by comparing phase
%           diff across cycles. This specifies the max gap in cycles for which a
%           comparsion can be made. Takes values from 1 plus where 1 implies only adjacent
%           cycles are analysed (i.e. 2 would allow a gap of 1 theta without spikes 
%           between the cyles with spikes [3].
%
% RETURNS
% spikePhase - [nSpikes x 1] for each spike its EEG phase (0 to 2pi) with spikes from
%           untrustworthy (low power) theta cycles set to nan
% intrinsicFreq - [nValidInterCyclePeriods x 1] intrinicFreq calculated for each semi-adjacent set of theta cycles that
%               contain spikes. Cyles that are more than maxElapsedCyles apart are
%               disgarded as are cycles and neighbours from low power cycles (defined by
%               validCycle).
%
% Original function


%%% Assign a phase to every spike %%%
instPhase250hz=mod(instPhase250hz,2*pi);            % Wrap phase (convert to phase between 0 and 2pi)
instPhase250hz(~ismember(cycle, validCycle))=nan;   % Set untrust worthy cycles to have phase nan
spikeSample250hz=ceil(spikeTimes.*sampleRate);      % Convert spiketime - measured in seconds by precise - to a specific 250Hz sample
spikePhase=instPhase250hz(spikeSample250hz);        % Spike phase - untrusworthy phases are set to nan

%%% Assign a cycle to every spike %%%
cycleForSpike=cycle(spikeSample250hz);          % Note: No attempt to deal with bursts that might cross between cycles
uniqueCycleSpike=unique(cycleForSpike);

%%% Get mean spike time for each cycle %%%
meanSample250hz=ones(size(uniqueCycleSpike))*nan;
for nn=1:length(uniqueCycleSpike)
%     meanPhase(nn)=circ_mean(spikePhase(cycleForSpike==uniqueCycleSpike(nn))); %The mean phase of spikes per cycle - nan for bad cyles
%     meanSample250hz(nn)=mean(spikeSample250hz(cycleForSpike==uniqueCycleSpike(nn))); % Mean spike time for cycle - units are EEG sample frequency counts.
    meanSample250hz(nn)=mean(spikeTimes(cycleForSpike==uniqueCycleSpike(nn))); % Mean spike time for cycle - units are seconds
end
% meanPhase=mod(meanPhase,2*pi); %Sort out fact that circ_mean gives neg values for >pi

%%% Get the elapsed time and cycles between the mean cycle spikes - these give the intrinsic frequency %%%
% elapsedTime=diff(meanSample250hz)./sampleRate;     % Time in secs between bursts
elapsedTime=diff(meanSample250hz);     % Time in secs between bursts (replaces line above when meanSample250hz calculated in sec, not samples)
apparentIntFreq=1./elapsedTime;                    % This is the 'absolute' intr freq: does not take into account how many cycles passed
elapsedCycles=diff(uniqueCycleSpike);              % Number of cycles between sampled cycles
intrinsicFreq=apparentIntFreq.*elapsedCycles;      % This divides the 'absolute' intr freq by N elapsed cycles. This gives us the intr freq in the 
                                                   % theta range. Note circularity here (how do we know cell is TM?). Maybe best to always look at adjacent cycles.

%%% Remove intrinsic frequencies where there are too many cycles between spikes, or there is an invalid cycle between spikes %%%
% Check for interposed invalid cycles %
interCyclesValid=zeros(size(intrinsicFreq));
for nn=1:length(interCyclesValid)
   range=uniqueCycleSpike(nn):uniqueCycleSpike(nn+1);
   interCyclesValid(nn)=all(ismember(range, validCycle));
end
% Check for n elapsed cycles %
intrinsicFreq=intrinsicFreq( elapsedCycles<=maxElapsedCycles & logical(interCyclesValid) );

% Put output in row vector format (so the histc always gives row vector output, even with 1 spike %
% if size(intrinsicFreq,1)size(intrinsicFreq,2) % If output is empty, want the opposite, i.e. size=(1,0)
if isempty(intrinsicFreq)
    intrinsicFreq=zeros(1,0);
    spikePhase=zeros(1,0);
else
    intrinsicFreq=intrinsicFreq';
    spikePhase=spikePhase';
end
% end
    



