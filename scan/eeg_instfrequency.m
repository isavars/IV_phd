function [instFreq,instPhase,instAmp] = eeg_instfrequency(eeg,passBand,varargin)
% Filters EEG & extracts instantaneous frequency, phase & power.
%   [instFreq,instPhase,instAmp] = eeg_instfrequency(data.trials.eeg,passBand)
%                                           OR
%   [instFreq,instPhase,instAmp] = eeg_instfrequency(eegData,passBand,sampRate)
%
% Note: phase of instPhase is monotonically increasing. Phase is not reset to 0 at 2pi but keeps
% increasing. Peak of EEG is defined as 0 phase (so peak is 0, 2pi, 4pi etc).

% Parse input %
if isstruct(eeg)
    sampleRate = eeg.sample_rate;
    eeg = ( double(eeg.eeg)./128 ) * eeg.scalemax;
elseif isnumeric(eeg)
    eeg = double(eeg);
    sampleRate = varargin{1};
end

% Filter EEG %
filteredEEG=eeg_filter(eeg,passBand,sampleRate);

% Calculate instantaneous parameters %
analyticEEG = hilbert(filteredEEG); %Hilbert transform itself - returns a complex var.

%Process analytic function to get phase, frequency and power (actually using true amp instead of power - power would just be that squared)
instPhase=angle(analyticEEG);       %Phase (radians) is the angle of the complex variable at each time point
instAmp=abs(analyticEEG);           %Modulus of analytic function i.e. instantaneous amplitude

%Calculate instantaneous freq - note this vector is length of EEG-1 to get
%round this add a nan to the end. Will result in very limited timeslipage
instFreq=diff(unwrap(instPhase))*(sampleRate/(2*pi));
instFreq(end+1)=nan;


