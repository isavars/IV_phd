function [rtn]=eeg_filter(eeg,passBand,varargin)
% Filter EEG signal.
%       filteredEEG=eeg_filter(eeg,passBand);            where eeg is a SCAn data.eeg struct.
%   	filteredEEG=eeg_filter(eegSamples,passBand,SR);  where eegSamples is eeg data, SR is sample rate.
%
% EEG is filtered with a 1-second long blackman window.
% Returned EEG values are in uV.

% Parse input %
if isstruct(eeg)
    eeg = ( double(eeg.eeg)./128 ) * eeg.scalemax;
    sampleRate = eeg.sample_rate;
elseif isnumeric(eeg)
    eeg = double(eeg);
    sampleRate = varargin{1};
end

% Filter EEG %
window = blackman(round(sampleRate)+1);
EEGFilter = fir1(round(sampleRate), passBand./(sampleRate/2) ,window);
rtn = filtfilt(EEGFilter, 1, eeg);