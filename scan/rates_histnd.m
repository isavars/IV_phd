function nd_array = rates_histnd(nd_data, nd_data_ranges, binsizes)
% HISTND Bins columnar N dimensional data into an Nd array. It uses a Matlab trick (sub2ind) to 
%   convert subscripts (sub) to indices (ind). These indices can then be binned very quickly 
%   using a 1D hist. The summed values are then assigned to Nd_array.
%
%   ND_ARRAY = RATES_HISTND(ND_DATA, ND_DATA_RANGES, BINSIZES)
%
%   Nd_data is a double array of size n_data_points x N (e.g. for spatial data, Nd_data = [x y dir])
%
%   NB: The lower limit of the range for each column of nd_data must be 0. If necessary, data 
%   should be transformed before calling this function (e.g. for xy spatial data, the transformation 
%   would be x_values = x_values - window_min_x; y_values = y_values - window_min_y). This shouldn't 
%   be needed for Axona output.
%
%   nd_data_ranges is a row vector of length N containing the maximum permitted range of the
%   data in the corresponding columns of Nd_data (e.g. for x-by-y-by-dir spatial data, 
%   nd_data_ranges = [(window_max_x - window_min_x) (window_max_y - window_min_y) (360 - 0)]
%
%   binsizes is a row vector of length N containing the size of the bins (in pixels) to be used 
%   to bin the corresponding columns of Nd_data (e.g. for spatial data, [8 8 6])

%   Written as a replacement for the Delaunay binning in mtint, Mike Anderson 30/8/09

[i,j] = size(nd_data);

% Check for correct input array sizes
% if size(nd_data_ranges,1) ~= 1 || size(nd_data_ranges,2) ~= j
%     error(' Ranges array is the wrong size.');
% end

% if size(binsizes,1) ~= 1 || size(binsizes,2) ~= j
%     error(' Binsizes array is the wrong size.');
% end

% Work out binned_array size
array_size = ceil(nd_data_ranges ./ binsizes);

% Check for empty input (no spikes) %
if isempty(nd_data)
    nd_array=zeros(array_size);
    return
end
    
% Form empty array
if j == 1
    nd_array = zeros(array_size,1);
else
    nd_array = zeros(array_size);
end

% Work out which bin each data point belongs to
binsizes = repmat(binsizes, [i 1]);
subs = ceil(double(nd_data) ./ binsizes);
if sum((subs == 0)~=0);
    disp('Zero indices in RATES_HISTND: check');
%     dbstop in rates_histnd at 57
    subs(subs == 0)=1;
end

% Convert subscripts to indices
if j == 1
    ind = subs;
elseif j ==2
    ind = sub2ind(array_size, subs(:,1), subs(:,2));
elseif j == 3;
    ind = sub2ind(array_size, subs(:,1), subs(:,2), subs(:,3));
else
    subs = num2cell(subs,1);
    ind = sub2ind(array_size, subs{:});
end

% Hist the indices across the total number of array bins
n_bins = prod(array_size);
binned_ind = hist(ind, [0.5:1:n_bins - 0.5]);









nd_array(1:n_bins) = binned_ind;


