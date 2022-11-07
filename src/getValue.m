%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trialInfo] = getValue(txt,keyStr)
    % Get a value from value-key cell array %
    ind = strcmp(keyStr,txt(:,1));
    if sum(ind) == 0
        trialInfo = [];
    else
        trialInfo = txt{ind,2};
    end
end