
function [rtn] = rates_mapmatch(dataVar, varargin)

% Find rate maps with specified paramaters in a data.rate_maps struct.
% If no match is found, rtn = [].
%
%       (1) index = rates_mapmatch(data, params);
%       (2) index = rates_mapmatch(data, 'paramname', value, .. etc );
%
% In the latter mode, all params not supplied are considered specified as
% default. 
% To check for a match regardless of a certain parameter, do
%
%       index = rates_mapmatch(data, .. , 'paramname', 'ignore');
%
% In calling mode 2 ONLY, 'crop', 'crop_to_filt' and 'abs_space' are
% ignored by default.
%
% First input argument can also be data.rate_maps.

% Sort out data input (for reasons I can't remember, can give either data
% or data.rate_maps.
if isfield(dataVar,'params')
    rateMaps = dataVar;
else
    rateMaps = dataVar.rate_maps;
end

% Check if any maps exist in input data struct %
if isempty(rateMaps)
    rtn = [];   return
end

% Parse input to make template params struct %
if isstruct(varargin{1})
    template = varargin{1};
else
    template = rates_params('default');
    template = rmfield(template, 'crop');
    template = rmfield(template, 'abs_space');
    template = rmfield(template, 'crop_to_filt');
    template = rmfield(template, 'filt_index');
    for ii=1:2:length(varargin)
        if strcmp(varargin{ii+1}, 'ignore')
            template = rmfield(template, varargin{ii});
        else
            template.(varargin{ii}) = varargin{ii+1};
        end
    end
end

% Check params in each rate_maps struct, check if all fields are the same %
fn = fieldnames(template);
rtn = [];
check=1;    % Now designed so that finding match is default (otherwise adding new fields to params struct caused havoc with this fucntion %
for ii=1:length(rateMaps)
    for jj=1:length(fn)
        f1 = template.(fn{jj});
        if isfield(rateMaps(ii).params,fn{jj})
            f2 = rateMaps(ii).params.(fn{jj});
            % Compare. This is a bit of a fuck-on, because I wasn't thinking of
            % doing this when I designed the params struct
            if ischar(f1)
                check = strcmp(f1, f2);
            elseif isempty(f1)
                check = isempty(f2);
            elseif all(size(f1)==size(f2))
                check = all(f1==f2);
            else
                check = 0;
            end
        end
        if ~check;   break;   end
    end
    if check;   rtn = ii;   break;   end
end

% if isempty(rtn);    
%     disp('RATES_MAPMATCH: no match found.');    
%     varargin
% end
            
        
    