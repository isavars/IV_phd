function [cycle, powerPerCycle, speedPerCycle] = eeg_powerspeedpercycle(eegPhase, eegAmp, speed)
% Identifies theta cycles and calculates power and speed per cycle. Used to identify cycles that 
% might contaminate analysis. Instant frequency and amplitude should previously have been calculated 
% for the EEG (filtered in theta band) using the Hilbert analytic signal. 
%
%   [cycle, powerPerCycle, speedPerCycle] = eeg_powerspeedpercycle(eegPhase, eegAmp, speed)
%
% ARGS
% eegPhase - instant phase for the EEG calculated using Hilbert. Must be 'unwrap'ed
% eegAmp - instant power for the EEG
% speed  - in cm/s, at pos sampling rate
%
% RETURNS
% cycle - [size(eegPhase)] vector indicating cycle that each sample point belongs to e.g.
%                   [0,0,1,1,1,1,1,2,2,2] etc
% powerPerCycle - [nx1] for each cycle the mean power (amplitude^2).
%                   NB have checked with Neil and this is power but not sure what units
%                   are. Need to check that with Ali.
% speedPerCycle - [nx1], mean speed per cycle.
%
% TW 04/03/10. Original function EEG_INSTANT_POWER_PER_CYCLE by CB, Jan '10.


% --- MAIN -------------------------------------------------------------------------
speedAtEEGSampRate=repmat(speed', 5, 1);                                    % col = 1 speed samp, row=1:5 reps of speed samp
speedAtEEGSampRate=reshape(speedAtEEGSampRate,numel(speedAtEEGSampRate),1); % Reshape to column vector

%NB. Analytic signal assigns 0rads to be the peak of the filtered signal. Therefore trough is pi.
phaseDiff = diff(mod(eegPhase, 2*pi)); %Change in phase between each pt - should be a saw tooth
phaseJump = find(phaseDiff<0); %Find point at which phase goes from 2pi to 0

nCycles=length(phaseJump)-1;
powerPerCycle=zeros(nCycles,1); %Pre allocate for speed
speedPerCycle=zeros(nCycles,1);
cycle=zeros(size(eegPhase)); %Vector same size as eegPhase
eegAmp=eegAmp.^2; %Square to turn amplitude to power

for nn = 1:nCycles
    speedPerCycle(nn) = mean(speedAtEEGSampRate( phaseJump(nn):phaseJump(nn+1) ));
    powerPerCycle(nn)=mean(eegAmp(phaseJump(nn):phaseJump(nn+1)));
    cycle(phaseJump(nn):phaseJump(nn+1))=nn;
end
