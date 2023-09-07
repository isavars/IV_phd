function [hMap] = gra_plotmap(map,env, varargin)

% Draw place map.
%       
%       [] = gra_plotmap(map);
%       [] = gra_plotmap(map, options);
%
% options.bm_thr  - Blank map threshold (default=1Hz)
%        .n_steps - Number of steps in map (default=5)
%        .scale_max - Max of step scale (default=peak rate)
%        .bw_colormap - greyscale map
%        .interp_uv - Flag: interpolate single unvisited bins? (default=0)
%        .text_pos - 'tl', 'none'
%        .handle - axes handle. (Default=[], draw into current axes).
%
% Can supply options as full struct, or as param/value pairs. In the latter
% case, unspecified options remain default.

% If no options supplied, generate defaults struct %
opt.bm_thr = 0;
opt.n_steps = 10;
opt.scale_max = [];
opt.bw_colormap = 0;
opt.interp_uv = 0;
opt.text_pos = 'tl';
opt.circ_mean = 0;  % Plot circular mean on polar HD plots.
opt.rotateDir=0;    % Rotate HD plots by this amount (in radians).
opt.handle = [];
opt.axes = [];
opt.parent = [];
if ~isempty(varargin) && ischar(varargin{1})
    for ii=1:2:length(varargin)
        opt.(varargin{ii}) = varargin{ii+1};
    end
elseif ~isempty(varargin) && isstruct(varargin{1})
    temp=varargin{1};  f=fieldnames(temp);
    for ii=1:length(f)
        opt.(f{ii}) = temp.(f{ii});
    end
end
if ~isempty(opt.parent);   opt.handle=opt.parent;   end
if ~isempty(opt.axes);   opt.handle=opt.axes;   end
map = double(map);

if isempty(opt.handle);   opt.handle = gca;   end
% IMAGE and PLOT wipe out axes tag, for some reason. This fix is necessary for SCAn map preview pane.
tagReset = get(opt.handle, 'tag');

%%%%%%%%%%%% Map plotting depends on place or direction map %%%%%%%%%%%%%%%
if size(map,2) > 1
    %%%%%%%%%%%%%%%%%%%%% Position Map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if isempty(map);   map = -1;   end   % Check for [] as input: draw a plain white image.
    map(isnan(map)) = -1;                % If background is NaN, convert to -1
    % Scale max? %
    if isempty(opt.scale_max)
        maxf = max(max(map));
    else
        maxf = opt.scale_max;
    end
    % Do interpolation, if wanted %
    if opt.interp_uv==1;   map = interp_uv(map);   end
    % Make indexed image: 2 modes. For 5-step maps, use a legacy method that
    % closely matches TINT. Otherwise use more general multistep routine.
    im = map;
    if opt.n_steps==5   
        im(map<0)       = 1;
        im(map>=0)      = 2;            
        if (maxf >= opt.bm_thr)                 % For blue map, stop here
            bdy = (round((maxf/5)*10)) / 10;
            im(map>bdy)     = 3;            % Values at boundary are placed in interval below, with 
            im(map>(bdy*2)) = 4;            % the exception of 0, which is rounded up to be in the lowest
            im(map>(bdy*3)) = 5;            % 'visited' interval.
            im(map>(bdy*4)) = 6;            %
        end
    else
        if maxf>opt.bm_thr
            im = im .* (opt.n_steps/maxf);
            im = ceil(im);
            im(im==0) = 1;    % Shift 0Hz rate to bin 1
        else
            im(:,:) = 1;        % This is to make a blue map.
        end
        im(map==-1) = 0;  % Unvis in bin 0
        im = im + 1;            % Indexed images must start at 1, not 0.
    end
    % Draw map %
    hMap = image(im, 'parent', opt.handle);
    if opt.bw_colormap
        g = linspace(0.7,0,opt.n_steps);
        colormap( cat(1, [1 1 1], repmat(g',1 ,3) ) );
    else
        colormap(opt.handle, gra_tintcolormap(opt.n_steps));
%         colormap(  cat(1,  [1 1 1;  jet(opt.n_steps)] )  );
    end
    axis(opt.handle, 'equal', 'tight', 'off');
    
elseif size(map,2) == 1
    %%%%%%%%%%%%%%%%%%%%% Direction Map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convert to x,y %
    bin_ang = (2*pi) / length(map);
    angles = (bin_ang/2 : bin_ang : (2*pi)-(bin_ang/2));
    angles=angles+opt.rotateDir;
    angles(angles>(2*pi))=angles(angles>(2*pi)) - 2*pi;
    angles(angles<=0)=angles(angles<=0) + 2*pi;
    [x, y] = pol2cart(angles', map);
    x(end+1) = x(1);    % Make line meet itself at end
    y(end+1) = y(1);    %
    % Plot %
    hMap = plot(opt.handle, x, y, 'k-');
    % Plot circular mean, if requested %
    if opt.circ_mean
        cm=circ_mean(angles',map);
        [cmX cmY]=pol2cart(cm,max(map));
        line('xdata',[0 cmX],'ydata',[0 cmY],'parent',opt.handle,'linewidth',2);
    end    
    % This is all formatting of axes %
    set(opt.handle, 'color', 'w');  % For SCAn preview pane.
    axis(opt.handle, 'tight');
    axlim = max(abs([get(opt.handle, 'xlim'), get(opt.handle, 'ylim')])); % Make axes square around biggest value
    axis(opt.handle,[-axlim axlim -axlim axlim]);                         %   ..
    line('xdata',0.95*[-axlim axlim],'ydata',[0 0],'parent',opt.handle);   line('xdata',[0 0],'ydata',0.95*[-axlim axlim],'parent',opt.handle); % centre-crossing axes
    axis(opt.handle, 'square', 'off', 'tight');
    
end
%%%%%%%%%%%%%%%%% Back to place/dir common code %%%%%%%%%%%%%%%%%%%%%%%%%%%

set(opt.handle, 'tag', tagReset);

% Max Rate Text %
rate_string = num2str( max(max(map)), '%4.1f') + " " + string(env) ; %'%4.1f
switch opt.text_pos     % Currently, you can have any position you want as long as it's top left.
    case 'tl'
        fs_val = 0.15;
        fu_val = 'normalized';
        pos_val = [0 1];
        hoz_val = 'left';
        ver_val = 'cap';
        text_color = [0 0 0];
        text('units', 'normalized', 'position', pos_val, 'HorizontalAlignment', hoz_val, 'string', rate_string, ...
            'FontUnits', fu_val, 'VerticalAlignment', ver_val, 'fontsize', fs_val, 'color', text_color,'parent',opt.handle);
    case 'bl'
        fs_val = 0.15;
        fu_val = 'normalized';
        pos_val = [1 1];
        hoz_val = 'left';
        ver_val = 'bottom';
        text_color = [0 0 0];
        text('position', pos_val, 'units', 'normalized', 'HorizontalAlignment', hoz_val, 'string', rate_string, ...
            'FontUnits', fu_val, 'VerticalAlignment', ver_val, 'fontsize', fs_val, 'color', text_color,'parent',opt.handle);
    case 'none'
        % No text
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rtn] = interp_uv(map)
% Interpolate values of single (8-surrounded) unvisited pixels %
% Columnise map (useful to pad first)  %
s = size(map);
pad = repmat(-1, s+2);
pad((2:(s(1)+1)), (2:(s(2)+1))) = map;
colmap = im2col(pad, [3 3], 'sliding');
% Find 8-surrounded unvisited bins %
ci = repmat(1, size(colmap));
ci(colmap==-1) = 0;
uvi = find(ci(5,:)==0);                     % unvisited index
s8i = find( sum( ci([1:4 6:9],:) )==8 );    % 8-surround index
u8i = intersect(uvi, s8i);
% Interpolate values, reshape to map %
colmap(5,u8i) = mean(colmap( [1:4 6:9], u8i ));
rtn = reshape(colmap(5,:), size(map));
