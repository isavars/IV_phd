function [mapsOut, varargout] = rates_main(data, varargin)
% Construct rate maps from raw data.
%
%       (1) R = rates_main(Data);
%       (2) R = rates_main(Data, 'paramName', paramValue, ... );
%       (3) R = rates_main(Data, Params);
%
% Where Data is a raw data strucure. Output R is a cell array of rate maps, {1:nTrial, 1:nCell}.
% Without any other inputs, makes a set of maps with default binning and smoothing parameters (1).
% To change these, you can either supply the parameter names and values as a comma-separated list
% of pairs (2). For the names and possible values of parameters, see HELP RATES_PARAMS.
% Or (3) you can directly pass the output of RATES_PARAMS (usually if you want to call the GUI associated 
% with that function).
%
% For some types of map there can be a more complex output:
%
%       [Rp, Rd]         = rates_main(data,Params);    (when param.alg is 'pxd' or 'dh')
%       [R, P, filtRads] = rates_main(data,Params);    (when using adaptive smoothing)


% Parse the input: is the secondary arg a struct of params or a name/value pair list? %
if isempty(varargin)
    params = rates_params('default');
elseif isstruct(varargin{1})
    params = varargin{1};
elseif ischar(varargin{1})
    params = rates_params(varargin{1:end}); 
end


% Check if we are making rate maps for all data or not %
if isempty(params.trial_index)
    trialsToMakeRateMaps = 1:length(data.trials);
else
    trialsToMakeRateMaps = params.trial_index;
end

% Preallocate ouput %
if strcmp(params.mode,'pos')
    mapsOut=cell(length(data.trials),1);
    mapsOut2=cell(length(data.trials),1);
else    
    mapsOut=cell(length(data.trials),length(data.trials(1).cells));
    if ~isempty(params.adaptive_smooth)
        mapsOut2=cell(length(data.trials),length(data.trials(1).cells));
    else
        mapsOut2=cell(length(data.trials),1);
    end
end

for ii=trialsToMakeRateMaps
    %%% Generate all filtered, binned position and spike maps %%
    % Get window (old SCAn compatibility) %
    if isfield(data.trials(ii), 'window_x')
        window = [ceil(data.trials(ii).window_y) ceil(data.trials(ii).window_x)];
    else
        window = [512 512];
    end
    % Correct wonky position values %
    data.trials(ii).dir(data.trials(ii).dir>359.9)=data.trials(ii).dir(data.trials(ii).dir>359.9) - 360; % Dirs should run 0-359. Sometimes 360s come out of POSTPROCESS_POS_DATA - dont have time to check why. (TW,21/05/08)
    data.trials(ii).dir=double(data.trials(ii).dir)+0.1;        % + 0.1 as CEIL should not create 0 indices.
    if min(data.trials(ii).x)==0;  data.trials(ii).x = data.trials(ii).x + 1;  end  % Once again, I've found some data with origin at 0,0. This shouldn't happen,
    if min(data.trials(ii).y)==0;  data.trials(ii).y = data.trials(ii).y + 1;  end  % but I don't have time to debug it properly. TW, 02/11/10.  
    data.trials(ii).x(data.trials(ii).x>window(2))=window(2);
    data.trials(ii).y(data.trials(ii).y>window(1))=window(1);
    % Crop (New crop: this won't be activated until clear of reproducibilty issues (25/09/09) %
%     [[data.trials(ii).x data.trials(ii).y]=rates_croppix(data.trials(ii).x,data.trials(ii).y,params.crop,50);
    % Filter positions %
    % Note that the filterPosData function now also removes x or y set to NaN (they, and corresponding spikes, are removed from map; function mod by TW, 09/09/2013).
    [filt_ind, filt_check] = filterPosData(data,ii,params);
    if sum(filt_ind) == 0 
        continue;           % If filter has left no pos samples, leave all in trial as [], go to next trial
    end
    % Set up bins depending on mode (saves time to have just one place bin for space='dir' etc) %
    if strcmp(params.alg,'pxd') || strcmp(params.alg,'dh')
        mapBins=[params.bin params.bin params.bin_dir];
    elseif strcmp(params.space,'place')
        mapBins=[params.bin params.bin 360];
    elseif strcmp(params.space,'dir')
        mapBins=[window params.bin_dir];
    end
    % Bin Position Map %
    pos_map=squeeze(rates_histnd([double(data.trials(ii).y(filt_ind)) double(data.trials(ii).x(filt_ind)) double(data.trials(ii).dir(filt_ind))],[window 360],mapBins));
    pos_map_temp=sum(pos_map,3);  %% NOTE: commenting this makes trouble for PxD?
    unvis=pos_map_temp==0;
    raw_maps{1}=pos_map./data.trials(ii).sample_rate;
    if strcmp(params.mode,'speed')
        % Speed maps %
        raw_maps{2} = rates_speedmap(data.trials(ii).y(filt_ind), data.trials(ii).x(filt_ind), data.trials(ii).dir(filt_ind), data.trials(ii).speed(filt_ind),window,mapBins);
    elseif ~strcmp(params.mode,'pos')
        %%% Spike maps %%%
        for jj=1:length(data.trials(ii).cells)
            % Convert times to positions samples %
            spk_samp=ceil(data.trials(ii).cells(jj).st .* data.trials(ii).sample_rate);
            % Select spikes by filter %
            spk_samp=spk_samp( ismember(spk_samp,find(filt_ind)) );
            % Spike map %
            raw_maps{jj+1}=squeeze(rates_histnd([double(data.trials(ii).y(spk_samp)) double(data.trials(ii).x(spk_samp)) double(data.trials(ii).dir(spk_samp))],[window 360],mapBins));
        end                                                                                          
    end
    
    % Filter binned maps for low occupancy, if requested %
    if ~isempty(params.min_dwell)
        lowOccInd = raw_maps{1}<=params.min_dwell;
        for jj=1:length(raw_maps)
            raw_maps{jj}(lowOccInd) = 0;
        end
        unvis = lowOccInd;
    end
    
    %%% Make Rate Maps %%
    if strcmp(params.alg,'pd')
        %%% (1) Regular rate Maps %%%
        if ~isempty(params.adaptive_smooth)
            if ~strcmp(params.mode,'rate')
                error('SCAn:rates_main','Adaptive smoothing is only applicable to mode=''rate''');
            end
            % Adaptive smoothing for place rate maps (Skaggs et al, 1996) %
            for jj=2:length(raw_maps)
                [mapsOut{ii,jj-1},~,mapsOut2{ii,jj-1},filtRads(ii,jj-1)]=rates_adaptivesmooth(raw_maps{1},raw_maps{jj},params.adaptive_smooth);
            end
        else
            % Regular p,d maps - rate, pos or spike %
            % Smooth %
            if strcmp(params.space,'place')           
                for jj=1:length(raw_maps)
                    raw_maps{jj}=smoothPlaceMap(raw_maps{jj},params.smooth,unvis);
                end    
            elseif strcmp(params.space,'dir')
                for jj=1:length(raw_maps)
                    raw_maps{jj}=smoothDirMap(squeeze(raw_maps{jj}),params.smooth_dir);
                end
            end
            % Assign to final output, making rate map if wanted %%
            if strcmp(params.mode,'pos')
                mapsOut(ii,1)=raw_maps(1);
            elseif strcmp(params.mode,'spike')
                mapsOut(ii,1:(length(raw_maps)-1))=raw_maps(2:end);
            elseif strcmp(params.mode,'speed')
                mapsOut(ii,1)=raw_maps(2);
            elseif strcmp(params.mode,'rate')
                for jj=2:length(raw_maps)
                    warning('off', 'MATLAB:divideByZero');
                    mapsOut{ii,jj-1}=raw_maps{jj}./raw_maps{1};
                    warning('on', 'MATLAB:divideByZero');
                end
                mapsOut2{ii,1} = raw_maps{1};
            end
        end
     
    elseif strcmp(params.alg,'pxd')
        %%% PxD %%%
        for jj=2:length(raw_maps)
            if sum(sum(sum(raw_maps{jj})))==0   % If no spikes, don't run algorithm.
                mapsOut{ii,jj-1}=zeros(size(raw_maps{jj},1),size(raw_maps{jj},2));   mapsOut{ii,jj-1}(unvis)=nan;
                mapsOut2{ii,jj-1}=zeros(size(raw_maps{jj},3),1);
            else
                [pMap, dMap, converged]=pxd_mtint(permute(raw_maps{jj},[3 1 2]),permute(raw_maps{1},[3 1 2]),30,0.00001,0.1); % pxd_mtint(spikes, times, max_iter, accuracy, tol)
                if converged
                    mapsOut{ii,jj-1}=smoothPlaceMap(pMap,params.smooth,unvis);
                    mapsOut2{ii,jj-1}=smoothDirMap(dMap,params.smooth_dir);
                end
                % NB. Non-converged maps are left as [].
            end
        end 
        
    elseif strcmp(params.alg,'dh')
        %%% Distributive hypothesis maps. rtn is dist hypo dir maps, rtn2 is dist hypo place maps. Output is already smoothed. %%%
        for jj=2:length(raw_maps)
            [mapsOut{ii,jj-1}, mapsOut2{ii,jj-1}]=rates_distributivehypothesis(raw_maps{1},raw_maps{jj},params.smooth_dir,params.smooth);
        end
        
    end
  
    % Crop (legacy crop - should be removed when croppix crop activated) %%
    if ~isempty(params.crop)
        if filt_check && length(params.crop)==2 && ~params.crop_to_filt
            % If pos is filtered, crop is to centre, and crop NOT to filt centres, make unfiltered pos maps to get crop coords.
            pos_map=rates_histnd([double(data.trials(ii).x) double(data.trials(ii).y) double(data.trials(ii).dir)],[window 360],[params.bin params.bin 360]);
            [~, crop_coords] = rates_crop(pos_map, ceil(params.crop/params.bin), 3);
        else
            crop_coords=ceil(params.crop/params.bin);
        end
        for jj=1:size(mapsOut,2)
            mapsOut{ii,jj} = rates_crop(mapsOut{ii,jj},crop_coords,3);
            if ~isempty(params.adaptive_smooth)
                mapsOut2{ii,jj} = rates_crop(mapsOut2{ii,jj},crop_coords,3);
            end
        end
    end
    
end

% Take only the rows (trials) required in output (if called with prms.trial_index~=[]) %
mapsOut = mapsOut(trialsToMakeRateMaps,:);
mapsOut2 = mapsOut2(trialsToMakeRateMaps,:);

% Assign varargout, if necessary %
if nargout==2
    varargout{1} = mapsOut2;
elseif nargout==3
    varargout{1} = mapsOut2;
    varargout{2} = filtRads;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sm_map]=smoothPlaceMap(map,k,unvis)
% Smooth a place map. map=map, k=boxcar kernel, unvis=index of unvisited bins.
if max(size(k))==1;  k=ones(k);  end   % Expand single parameter to flat k-by-k square
map(unvis)=0;
visTemplate=ones(size(map));
visTemplate(unvis)=0;
filtMap=imfilter(map,k);
filtVis=imfilter(visTemplate,k);
warning('off', 'MATLAB:divideByZero');
sm_map=filtMap./filtVis;
warning('on', 'MATLAB:divideByZero');
sm_map(unvis)=nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sm_map]=smoothDirMap(map,window)
% Smooth a dir map.
if window==1; sm_map=map; return; end
p = (window-1)/2;                                               % Pad for circular smooth
pad_map = [map(end-p+1:end); map; map(1:p)];                    %  ..
sm_map = mean( im2col(pad_map, [window 1], 'sliding') )';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [filt_ind, filt_active] = filterPosData(data, tr_ind, params)
% Produces an index of the position samples allowed
% through by filters.
%
%       [filt_ind, filt_active] = filterPosData(data, trial_ind, params);
%
% If no filters are active, or no points are excluded, 
% filt_ind will be 1:length(data.trials(trial_ind).x).
% filt_active is a flag for if the data has been filtered.

filt_ind = true(size(data.trials(tr_ind).x));
filt_active = 0;

% Direct index filter %
if isfield(params, 'filt_index') && ~isempty(params.filt_index)
    if iscell(params.filt_index)
        filt_ind = filt_ind & params.filt_index{tr_ind};
    elseif islogical(params.filt_index) && length(params.filt_index) > length(data.trials(tr_ind).x)
        filt_ind = filt_ind & params.filt_index( 1:length(data.trials(tr_ind).x) );
    elseif islogical(params.filt_index) && length(params.filt_index) == length(data.trials(tr_ind).x)
        filt_ind = filt_ind & params.filt_index;
    else
        error('Direct index filter parameters in incorrect format. Should be 1xnPosSample or {nTrial}x1xnPosSample.'); 
    end
    filt_active = 1;
end

% Filter 1st or 2nd half of the trial %
if params.filt_1stHalf
    tempInd=false(size(data.trials(tr_ind).x));
    tempInd(   1:round(data.trials(tr_ind).dur/2 * data.trials(tr_ind).sample_rate)   ) = true;
    filt_ind = filt_ind & tempInd;
end
if params.filt_2ndHalf
    tempInd=false(size(data.trials(tr_ind).x));
    tempInd(   round(data.trials(tr_ind).dur/2 * data.trials(tr_ind).sample_rate) : length(data.trials(tr_ind).x) ) = true;
    filt_ind = filt_ind & tempInd;
end

% Filter by time %      
if ~isempty(params.filt_time)
    filt_active = 1;
    % NOTE can give params.filt_times as n_trialx2, one param for each trial
    if all(size(params.filt_time) == [1 2])
        time_param = params.filt_time;
    elseif all(size(params.filt_time) == [length(data.trials), 2])
        time_param = params.filt_time(tr_ind,:);
    else
        error('Time filter parameters in incorrect format. Should be 1x2 or n_trialx2.'); 
    end    
    samp_rate = data.trials(tr_ind).sample_rate; 
    filt_index_start = round((time_param(1) * samp_rate) + 1);
    filt_index_end = round(time_param(2) * samp_rate);
    tempInd=false(size(data.trials(tr_ind).x));
    tempInd(filt_index_start:filt_index_end) = true;
    if length(tempInd) ~= length(filt_ind)
        tempInd = tempInd(1:length(filt_ind));  % If the time filter parameters run over the actual length of the trial, restrict the length of the time filter ind to the trial length. 
    end
    filt_ind = filt_ind & tempInd;
end

% Filter by speed %
if ~isempty(params.filt_speed)
    filt_active = 1;
    tempInd = data.trials(tr_ind).speed > params.filt_speed(1) & data.trials(tr_ind).speed <= params.filt_speed(2);
    filt_ind = filt_ind & tempInd;
end

% Filter by direction %
if ~isempty(params.filt_dir)
    filt_active = 1;
    temp1 = data.trials(tr_ind).dir >= params.filt_dir(1);
    temp2 = data.trials(tr_ind).dir <= params.filt_dir(2);
    % Deal with circularity here. Always select anti-clock from value 1 to value 2
    if params.filt_dir(1) < params.filt_dir(2)
        dir_filt = temp1 & temp2;        
    elseif params.filt_dir(1) > params.filt_dir(2)
        dir_filt = temp1 | temp2;
    end   
    filt_ind = filt_ind & dir_filt;
end

% Filter by position in X and/or Y %
if ~isempty(params.filt_x)
    filt_active = 1;
    tempInd = data.trials(tr_ind).x > params.filt_x(1) & data.trials(tr_ind).x <= params.filt_x(2);
    filt_ind = filt_ind & tempInd;
end
if ~isempty(params.filt_y)
    filt_active = 1;
    tempInd = data.trials(tr_ind).y > params.filt_y(1) & data.trials(tr_ind).y <= params.filt_y(2);
    filt_ind = filt_ind & tempInd;
end

% Filter by position = NaN (new code, TW, 09/09/13) %
if sum(isnan(data.trials(tr_ind).x)) > 0   ||   sum(isnan(data.trials(tr_ind).y)) > 0 
    filt_active = 1;
    tempInd = ~isnan(data.trials(tr_ind).x)   &   ~isnan(data.trials(tr_ind).y);
    filt_ind = filt_ind & tempInd;
end

















