function [rtn] = rates_main_old_NIU(data, params)

% NOTE, 01/10/09. This version was superceded by PxD capable version.

% Construct rate maps from raw data.
%
%       M = rates_main(data, Params);
%
% Where data is a data strucure, Params is structure giving rate map parameters.
%
% Output M is a cell array of rate maps, {1:nTrial, 1:nCell}.

for ii=1:length(data.trials)
%     % Half Maps Hack %
%     params.filt.time = [round(data.trials(i).dur/2) data.trials(i).dur];
%     data.rate_maps(map_ref).name = 'rate-place: 2nd Half, no crop';
    
    % Get window (old SCAn compatibility) %
    if isfield(data.trials(ii), 'window_x')
        window = [data.trials(ii).window_x data.trials(ii).window_y];
    else
        window = [512 512];
    end
   
	samp_rate = data.trials(ii).sample_rate;
    x_vector = double(data.trials(ii).x);
    y_vector = double(data.trials(ii).y);

    %% Filters %%
    [filt_ind, filt_check] = filterPosData(data, ii, params);
    if isempty(filt_ind)    
        continue;           % If filter has left no pos samples, leave all in trial as [], go to next trial
    end
    x_vector = x_vector(filt_ind);
    y_vector = y_vector(filt_ind);
    
    %%% Position Map %%%
    if strcmp(params.space, 'place')
        [pos_map, x_bin, y_bin] = rates_binpos(x_vector, y_vector, params.bin, window);
    elseif strcmp(params.space, 'dir')
        [pos_map, dir_bin] = rates_binpos(double(data.trials(ii).dir(filt_ind)), params.bin);
    end
    pos_map = pos_map ./ samp_rate;
    
    % Crop %
    if ~isempty(params.crop)
        cropped_size = ceil(params.crop/params.bin);
        if filt_check && length(cropped_size)==2 && ~params.crop_to_filt
            % If the positions have been filtered, crop is to centre, and crop to filt centres
            % not asked for make whole trial positions as the basis for cropping.
            [whole_pos_map, junk, junk2] = ...
                rates_binpos(double(data.trials(ii).x), double(data.trials(ii).y), params.bin, window);   
            [junk, rtn_crop_coords] = rates_crop(whole_pos_map, cropped_size, 3);
            [pos_map, junk] = rates_crop(pos_map, rtn_crop_coords, 3);
        else
            [pos_map, rtn_crop_coords] = rates_crop(pos_map, cropped_size, 3);
        end
    end   
    unvis = find(pos_map==0);
    
    % Smooth %
    if params.smooth > 1  %% Need to put this check in because colfilt gives an error for neighbourhood (smoothing) = 1.
        if strcmp(params.space, 'place')
            pos_map(unvis) = NaN;
            pos_map = colfilt(pos_map,[params.smooth,params.smooth],'sliding',@nanmean);
        else
            p = (params.smooth - 1) / 2;                              % Pad for circular smooth
            pad_map = [pos_map(end-p+1:end); pos_map; pos_map(1:p)];  %  ..
            pos_map = mean(im2col(pad_map, [params.smooth 1], 'sliding'))';
        end
    end
    pos_map(unvis) = NaN;
    
    if strcmp(params.mode, 'pos');
        % Output pos map only %
        rtn{ii,1} = single(pos_map);
    else
        for jj=1:length(data.trials(ii).cells)
            %%% Spike maps %%%
            % Convert spike times into spike count vector ([1:n_pos_samp], index into position samples).
            % Doing this, then going back individual spikes to make the spike map (below)
            % is inefficient, but makes selecting spikes by filter much easier.
            spk_samp = data.trials(ii).cells(jj).st .* samp_rate;             
            spk_samp = floor(spk_samp) + 1;                         
            spk_vect = repmat(0, [1 length(data.trials(ii).x)]);
            for kk=1:length(spk_samp);
                spk_vect(spk_samp(kk)) = spk_vect(spk_samp(kk)) + 1;       
            end                                                           
            
            % Select spikes by filter %
            spk_vect = spk_vect(filt_ind);

            % Spike vector -> spike map %
            if strcmp(params.space, 'place')
                mapsize = ceil(fliplr(window)/ params.bin);  % Note x,y to r,c conversion
                spk_map = repmat(0, mapsize);
                spk_fires = find(spk_vect);
                for kk=1:length(spk_fires)                  
                     curr_pos_col = x_bin(spk_fires(kk));
                     curr_pos_row = y_bin(spk_fires(kk));
                     spk_map(curr_pos_row, curr_pos_col) = spk_map(curr_pos_row, curr_pos_col) + spk_vect(spk_fires(kk));
                end
            elseif strcmp(params.space, 'dir')
                spk_fires = find(spk_vect);
                spk_map = repmat(0, length(pos_map), 1);
                for kk=1:length(spk_fires)                  
                     curr_pos = dir_bin(spk_fires(kk));
                     spk_map(curr_pos) = spk_map(curr_pos) + spk_vect(spk_fires(kk));
                end
            end
            
            % Crop, Smooth, Rate %
            if ~isempty(params.crop);
                spk_map = rates_crop(spk_map, rtn_crop_coords);
            end
	        
            if params.smooth > 1
                if strcmp(params.space, 'place')
                    spk_map(unvis) = NaN;
                    spk_map = colfilt(spk_map,[params.smooth,params.smooth],'sliding',@nanmean);
                else
                    p = (params.smooth - 1) / 2;                              % Pad for circular smooth
                    pad_map = [spk_map(end-p+1:end); spk_map; spk_map(1:p)];  %  ..
                    spk_map = mean(im2col(pad_map, [params.smooth 1], 'sliding'))';
                end
            end
            
            % Spike map or rate map? %
            if strcmp(params.mode, 'spike')
                output_map = spk_map;
            else
                warning('off', 'MATLAB:divideByZero');    % Background 0's - go to NaN anyway.
                output_map = spk_map ./ pos_map;
                warning('on', 'MATLAB:divideByZero');
            end  
            output_map(unvis) = NaN;
            
            % assign map in output %
            rtn{ii,jj} = single(output_map);
 
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

filt_ind = 1:length(data.trials(tr_ind).x);
filt_active = 0;

% Direct index filter %
if isfield(params, 'filt_index') && ~isempty(params.filt_index)
    if size(params.filt_index,1)==1
        params.filt_index=repmat(params.filt_index,length(data.trials),1);
    end
    for ii=1:length(data.trials)
        filt_ind = intersect(params.filt_index(ii,:),filt_ind);
    end
    filt_active = 1;
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
    filt_ind = intersect(filt_ind, (filt_index_start:filt_index_end));
end

% Filter by speed %
if ~isempty(params.filt_speed)
    filt_active = 1;
    spd_cms = params.filt_speed * (data.trials(tr_ind).ppm / 100);    % To convert the supplied parameter, in cm, to pixels
    faster_min = find( data.trials(tr_ind).speed > spd_cms(1) );
    slower_max = find( data.trials(tr_ind).speed < spd_cms(2) );
    spd_filt = intersect(faster_min, slower_max);
    filt_ind = intersect(filt_ind, spd_filt);
end

% Filter by direction %
if ~isempty(params.filt_dir)
    filt_active = 1;
    temp1 = find( data.trials(tr_ind).dir >= params.filt_dir(1) );
    temp2 = find( data.trials(tr_ind).dir <= params.filt_dir(2) );
    % Deal with circularity here. Always select anti-clock from value 1 to value 2
    if params.filt_dir(1) < params.filt_dir(2)
        dir_filt = intersect(temp1, temp2);        
    elseif params.filt_dir(1) > params.filt_dir(2)
        dir_filt = union(temp1, temp2);
    end   
    filt_ind = intersect(filt_ind, dir_filt);
end

















