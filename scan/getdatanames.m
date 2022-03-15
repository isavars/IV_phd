function [rtn] = getdatanames(data, type, option, varargin)

% Get a list of formatted cell or trial names. (In cell array of strings).
%
%       [rtn] = getdatanames(data, 'type', 'option');
%       [rtn] = getdatanames(data, 'type', 'option', index);
%
% 'type' - 'cell' or 'trial'
% 'option' - 'long' or 'short'
% For cells, 'Cell 12 Tet 4' (long) or 'c12 t4' (short).
% For trials, '081207a' (long) or 'a' (short).
%
% index will return only a selection of the names, based on the index.

switch type
    case 'trial'
    	for ii=1:length(data.trials) 
        	if strcmp(option,'short')
                rtn{ii} = data.trials(ii).trialname(end);
            else
                rtn{ii} = data.trials(ii).trialname;
            end
        end
	case 'cell'
        for ii=1:length(data.trials(1).cells)
            if strcmp(option,'long')
                rtn{ii} = ['Cell ', num2str(data.trials(1).cells(ii).cellnum), ' Tet ', num2str(data.trials(1).cells(ii).tet)];
            else
                rtn{ii} = ['c', num2str(data.trials(1).cells(ii).cellnum), ' t', num2str(data.trials(1).cells(ii).tet)];
            end
        end
end

if ~isempty(varargin)
    ind = varargin{1};
    rtn = rtn(ind);
end
        