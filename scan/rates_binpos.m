function [varargout] = rates_binpos(varargin)

% Bin position-by-sample vectors, both place and direction data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For place data:
%
%     [pos_map, x_binned, y_binned] = rates_binpos(x, y, binsize, [window_x window_y]);
%
% x and y are the vectors of position samples.
% binsize is the side length, in pixels, of a bin.
% [window_x window_y] is the size of the input image, in pixels. (Need 
% this so map is correctly 'framed' in window in output).
%
% pos_map is a 2-D array showing positions visited, but careful, units are 
% absolute numbers of dwell samples. Need to divide by position sample rate
% to give pos map in seconds units.
%
% x_binned and y_binned are the binned x,y vectors, still in vector form. 
% ie. if x(n) = a, x_binned(n) = ceil(a/binsize). (Use these to make binned spike maps).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For direction data:
%
%     [pos_vect, dir_bin] = rates_binpos(dir, binsize);
%
% Binsize should be specified in degrees, and there should go into 360. 
% Bins run 0 -> binsize-1, binsize -> binsize*2-1, etc. (As TINT convention 
% has directions going 0-359).
%
% pos_vect is a vector containing the number of position dwell
% samples at each heading ([0:binsize:359]).
%
% dir_bin is the input directional position vector, binned.
% Headings in the range 0 - (binsize-1) have value 1, 
% binsize-(binsize*2-1) have value 2, etc.

if length(varargin) == 4
    %%%%%%%%%%%%%%%%%%%%%%%% Place Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = varargin{1};
    y = varargin{2};
    binsize = varargin{3};
    window = varargin{4};
    
    % Clip data out of window %
    x(x>window(1)) = window(1);
    y(y>window(2)) = window(2);
    
    %%%  Bin Image %%%
    xb = ceil(x/binsize);             % Bin the vector by doing (/8 and round) to x and y vectors. 
    yb = ceil(y/binsize);             % Use ceil, so bins go 1-8, 9-16 etc.
    winBin = ceil(fliplr(window) ./ binsize);       % Binned map size (round up any 'half bins' at the edge to whole bins). FLIPLR: x,y to r,c.

    %%% Prepare xy for resampling to matrix image %%%
    ind = sub2ind(winBin, yb, xb);             % x and y to linear index. NOTE: x and y switched, to convert real-world x,y to matlab rows,columns
    ind = sort(ind);                           % Sort - puts all values for each position together, in linear index order.                              
    binstep = find(diff(ind));                 % Find where index changes, ie boundaries for values between bins.
    binstep = [0; binstep; length(ind)];       % Add bounding for first and last pixel.

    %%% Create binned image %%%
    binstep_length = length(binstep);                   % Take this out of loop for speed increase (over-optimistic?).
    binvect = zeros(winBin(1)*winBin(2),1);             % Pre-allocate binned map vector
    n_corr = 1;                                         % n_corr will define where to look in binstep (and therefore ind). It will only increase in step with n when ind contains entries relevant to binvect(n).
    for ii=1:(winBin(1)*winBin(2))       % for every bin in output ...      
        if n_corr==binstep_length                       % Check if we have reached the end of binstep (and therefore ind) ...
            break                                       % ... if so, leave the rest of binvect as 0
        end                                             % (This bit deals with empty bottom right corners) 

        if ind(binstep(n_corr)+1)==ii                    % Does the linear index value in ind (binstep(n_corr)+1) correspond with that in binvect (n) ?        
            binvect(ii) = binstep(n_corr+1) - binstep(n_corr);          % the difference between adjacent binstep values give the number of times the animal has dwelt in this bin
            n_corr = n_corr+1;                                          % .. and increase n_corr
        end
    end
    binmap = reshape(binvect,winBin(1),winBin(2));      % Convert binned vector to 2D matrix.
    % Output %
    varargout{1} = binmap;
    varargout{2} = xb;
    varargout{3} = yb;
       
elseif length(varargin) == 2
    %%%%%%%%%%%%%%%%%%%%% Direction Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Input %
    dir = varargin{1};
    binsize = varargin{2};
    if rem(360,binsize)~=0
        error('Binsize doesn''t go into 360 degrees');
    end
    % Bin directions %
    dir = double(dir);
    dir(dir==360)=0;  % Dirs should run 0-359. Sometimes 360s come out of POSTPROCESS_POS_DATA - dont have time to check why. (TW,21/05/08)
    dir_bin = ceil((dir+1)./binsize);    % + 1 as dir is 0-359, not 1-360.
    % Convert to position vector (summed dwell times) %
    sort_bin = sort(dir_bin);                     % Sort - puts all values for each position together, in linear index order.                              
    binstep = find(diff(sort_bin));               % Find where index changes, ie boundaries for values between bins.
    binstep = [0; binstep; length(sort_bin)];     % Add bounding for first and last pixel.
    pos_vect = diff(binstep);
    % Assign Output %
    varargout{1} = pos_vect;
    varargout{2} = dir_bin;
    
end
    
    
    
    
    
    
    
