function [eegInd]=eeg_posind2eegind(posInd)
% Convert a (numerical) index into positions into the matching index into EEG data.
% Asssumes EEG sample rate 5x Pos sample rate. (i.e. no good for egf data).
%       [eegInd]=eeg_posind2eegind(posInd)
if size(posInd,1)==1
    posInd=posInd';
end
a = repmat(posInd.*5, 1, 5);
b = repmat(0:4, length(posInd), 1);
eegInd= a-b;
eegInd = sort(reshape(eegInd, [], 1));