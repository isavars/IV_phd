
function [cropped_map,crop_out] = rates_crop(map,coords,e)

% Find centre and crop
% 
%           M = rates_crop(map, coords, e);
%      [M, C] = rates_crop(map, coords, e);
%
% coords format (3 options):
%
%  [A B]     - Crop a box, AxB (width x height), around the environment centre.
%  [x y A B] - Crop absolute co-ordinates, defined [top_left_x, top_left_y, width, height].
%  'tight'   - Crop the map tight to the extent of the visited bins.
%
% e = Neil's edge trimming variable. Set 0 to ?. Not relevant for absolute coords, but need to give [] as placeholder.
% M = cropped map matrix.
% C = return of crop co-ordinates used, in format [top_left_x, top_left_y, width, height].

% NOTE. In and out co-ord args are in format [Width, Height] (or [x, y, width, height]).
%       However, internally (from line 33 to line 115) this program uses format [rows, columns].

% Check coords format. (Also switch x,y to row, column) %
if ischar(coords)
    find_edges = 1;
    find_cen = 0;
elseif all(size(coords) == [1 2])
    find_edges = 1;
    find_cen = 1;
    coords = fliplr(coords);
elseif all(size(coords) == [1 4])
    find_edges = 0;
    coords(1:2) = fliplr(coords(1:2));
    coords(3:4) = fliplr(coords(3:4));
else
    error('Error in RATES_CROP. Crop co-ords in wrong format.');
end

[size_row, size_col] = size(map);

if min([size_row,size_col])==0
    cropped_map=map;
    crop_out=[];
    return
end

if find_edges == 1
    % Find where edges of visited enviroment are %
    if sum(sum(isnan(map)))>0                   % If background is NaN.
        a = find(~isnan(map));                  % find visited bins (for rows)
		ainv = find(~isnan(map'));              % find visited bins (for columns)
        bkgnd=nan;
	else                                       % if 0..
        a = find(map);                         % same again
		ainv = find(map');                     %  ..
        bkgnd=0;
    end
	if (e+1) > length(a)
        e = length(a) - 1;  % In case of very small visited enviroments.
	end
	col = [ceil(a(e+1)./size_row), ceil(a(end-e)./size_row)];               % find edge defining columns (NOTE: 'size_row' denotes the number of rows, i.e. the length of each column).
	row = [ceil(ainv(e+1)./size_col), ceil(ainv(end-e)./size_col)];         % find edge defining rows    (NOTE: 'size_col' denotes the number of columns, i.e. the length of each row).
	
    if find_cen == 1
        % Find centre %%
        cen = floor([mean(row), mean(col)]);    % round down, then do (d-1) crop in the negative-going direction, below.
        % Make Crop Indices %%
        offset_pos_row = floor(coords(1)/2);
        offset_pos_col = floor(coords(2)/2);
        offset_neg_row = ceil(coords(1)/2) - 1;
        offset_neg_col = ceil(coords(2)/2) - 1;

        crop(1) = cen(1) - offset_neg_row;      % Top row 
        crop(2) = cen(2) - offset_neg_col;      % Left col 
        crop(3) = cen(1) + offset_pos_row;      % Bottom row
        crop(4) = cen(2) + offset_pos_col;      % Right col       
    elseif find_cen == 0
        % Crop tight to edges %%
        crop(1) = row(1);      % Top row 
        crop(2) = col(1);      % Left col
        crop(3) = row(2);      % Bottom row
        crop(4) = col(2);      % Right col      
    end
    
elseif find_edges == 0
    % Predefined co-ords %%
    crop(1) = coords(1);
    crop(2) = coords(2);
    crop(3) = crop(1) + coords(3) - 1;
	crop(4) = crop(2) + coords(4) - 1;
end

% In case the crop co-ordinates go 'off the map', pad map with background %
pad = 0;
if any(crop < 1)
    neg_bound = min(crop);
    pad = 1;
else
    neg_bound = 1;
end
if (crop(3) > size_row) || (crop(4) > size_col)
    pos_bound = max(crop);
    pad = 1;
else
    pos_bound = max(crop(3:4));
end
if pad == 1
    % pad map %
    padded = repmat(nan, (pos_bound-neg_bound+1), (pos_bound-neg_bound+1));
    pad_tl_corner = (neg_bound * -1) + 1;
    padded( (1:size_row)+pad_tl_corner, (1:size_col)+pad_tl_corner ) = map;
    map = padded;
    % adjust crop %
    pad_crop = crop + (1-neg_bound);
else
    pad_crop = crop;
end

% Crop %
cropped_map = map(pad_crop(1):pad_crop(3), pad_crop(2):pad_crop(4));           % Crop around centre.

% Make crop coords return arg. Note switch back from row,col to x,y %
crop_out(1:2) = fliplr(crop(1:2));
crop_out(3) = crop(4) - crop(2) + 1;
crop_out(4) = crop(3) - crop(1) + 1;

