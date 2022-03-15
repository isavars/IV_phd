function [transMapsTransformed] = rates_transform(transMaps,shapeCodeTrans,targetMap,shapeCodeTarget,varargin)
% Produce topologically transformed maps. Registering size is optional, and can rotate to an extent (see below).
% Works for squares, circles, and octagons, NOT for square/plus data.
%
%    [mapsTransformed] = rates_transform(transMaps,shapeCodeTrans,targetMap,shapeCodeTarget, .. options .. )
%
% transMaps is a {1,nCell} cell array of all rate maps from the trial to be transformed.
% targetMap is a rate map from a trial of the target shape.
%
% shapeCodeX is the shape definer, 8,7,6,5,4,1 (Square, 1:7 Oct, 6:2 Oct, 5:3 Oct, 4:4 Oct, Circle)
%
% Options: 'reg_size', {1, 0}   - Register environment size    (Default = 1)
%          'rotate', D          - Rotate shape 2 by anti-clock D degrees. (*[tr_x, tr_y]* rotated 
%                                 clockwise; after INTERP2-ing, *shape 2 rates* rotated anticlockwise).
%          'force_scale', SF    - forces a particular scaling factor, SF, onto the transformed maps.
%
% IMPORTANT NOTE ON ROTATIONS. Rotations are done after the transformation - they won't
% transform a rotated shape. Rotations will only give meaningful results, therefore, if 
% the rotation specified is a rotational symmetry of the shape being transformed (shape 2).
%
% Transformation is based on 1) registering shapes on the basis of diameter and 
% centroid position, 2) radial transformation of bin position, based on relative 
% radii size of shape perimeters at any given angle, 3) nearest-neighbour 2-D 
% interpolation from original rate bin values to transformed position ones.


% Input options %
prms.reg_size = 1;
prms.rot = 0;
prms.force_scale = 0;
for ii=1:2:length(varargin)
    prms.(varargin{ii}) = varargin{ii+1};
end
if ischar(shapeCodeTarget); shapeCodeTarget=str2double(shapeCodeTarget); end  % If shapes have been passed from SCAn .user variables
if ischar(shapeCodeTrans); shapeCodeTrans=str2double(shapeCodeTrans); end     % without converison from char.

[xx, yy] = meshgrid(1:size(targetMap,1), 1:size(targetMap,2));

% Find environment size and position %
E1=trans_findenv(targetMap,shapeCodeTarget);   % TRANS_FINDENV is a subfunction.
E2=trans_findenv(transMaps{1},shapeCodeTrans);

% Get x,y coords of bins in Env, and convert to polar coords, (0,0) in centre of environment %
x = xx(E1.in) - E1.centre(1);
y = yy(E1.in) - E1.centre(2);
[th, r] = cart2pol(x,y);

% Transform shape radially %
r_tr = trans_radcoords(th, r, shapeCodeTarget, shapeCodeTrans);    % TRANS_RADCOORDS is a subfunction

% Register size %
if prms.force_scale
    r_tr = r_tr ./ prms.force_scale;    % Note that this is the inverse of the given scaling factor (as r_tr is the transform of the target to the transformed).
elseif prms.reg_size == 1
	% first, need to check what size relationship should be, given shapes 
	shape_defs = [shapeCodeTarget, shapeCodeTrans];
    rad_temp = [0 0];
	for i=1:2
        switch shape_defs(i)
        case 1                  % Circle
            rad_temp(i) = 1;
        case 8                  % Square
            rad_temp(i) = pi/4;
        otherwise               % Octagons
            n = shape_defs(i);
            rad_temp(i) = (((1+sqrt(2))*pi)./8) *  ((n/2) + sqrt(((8-n)^2)/2)) * (1/(2+sqrt(8)));    % see TRANS_RADCOORDS for details.
        end
	end
	shape_ratio = rad_temp(1) / rad_temp(2);  % This is what the radius ratio should be.
	actual_ratio = ((E1.box(3)+E1.box(4))/4) / ((E2.box(3)+E2.box(4))/4);  % This is what it is.
	r_tr = r_tr .* (shape_ratio/actual_ratio);
end

% Rotate %
if prms.rot
    prms.rot = (prms.rot/360) * 2 * pi;
    th = th + prms.rot;
end

% Back to xy coords, register centroid %
[xTr,yTr] = pol2cart(th, r_tr);
xTr = xTr + E1.centre(1) + (E2.centre(1) - E1.centre(1));  % Move the shape1 xy coords to be 
yTr = yTr + E1.centre(2) + (E2.centre(2) - E1.centre(2));  % centred at the shape2 centre.

% Transform Rate Maps, by interpolating values for transformed x,y coords from to-be-transformed rate maps %
for ii=1:length(transMaps)
    % First do 2-D interpolation %
    transTemp = interp2(xx, yy, transMaps{ii}, xTr, yTr, 'nearest')';
    % Now reconstitute 'map', i.e. 2-D array corresponding to place %
    mapTemp = nan(size(targetMap));
    mapTemp(E1.in) = transTemp;
    transMapsTransformed{ii}=mapTemp; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Env] = trans_findenv(times, shape)
% Find edges of visited environment.
%
%       Env = trans_findenv(times, shape);
% 
% Env.centre = [x y]
% Env.box = [x y width height]
% Env.in = list of indices inside the occupied area
% 
% uses image processing toolbox

occ = zeros(size(times)); 
occ(~isnan(times)) = 1;

occ = bwmorph(occ,'clean');
occ = bwmorph(occ,'spur', Inf);

[y x] = find(occ > 0);
[xx, yy] = meshgrid(1:size(times,1), 1:size(times,2));

k = convhull(x, y);  
filled_occ = inpolygon( yy, xx, y(k), x(k));           
in = find(filled_occ>0);            % works for convex environemnts
filled_occ(in) = 1;                 % points on the boundary are given as 0.5 - include in list of indeces
stats = regionprops(double(filled_occ), 'Centroid', 'BoundingBox');      % BoundingBox & Centroid OK if concave.

if(shape > 17)        % concave envs! - plus mazes. 
    stats = regionprops(occ, 'Centroid', 'BoundingBox', 'FilledImage');
    % remove points in filled_occ that are unfilled in equiv FilledImage bin (NB FilledImage is cropped).
    xstartpad = min(xx(in))-1;
    ystartpad = min(yy(in))-1;
    xendpad = size(xx, 2) - max(xx(in));
    yendpad = size(xx, 1) - max(yy(in));
    % pad out FilledImage to match xx, yy, filled_occ, in etc
    temp = [ zeros(size(stats.FilledImage, 1), xstartpad) stats.FilledImage zeros(size(stats.FilledImage, 1), xendpad) ];
    temp = [ zeros( ystartpad, size(xx, 2) ); temp; zeros( yendpad, size(xx, 1) ) ];
    filled_occ( temp < 1 ) = 0;
    in = find( filled_occ > 0 );
end

% TW, 12/12/06. Output repackaged for my convenience %
Env.centre = stats.Centroid;
Env.box = stats.BoundingBox;
Env.in = in;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rT] = trans_radcoords(th, r, shape_in, shape_out)
% Produce transformed radial coords (square, circle and octagon) 
% on basis of radial distance from centre.
%
%   [rT] = trans_radcoords(Th, R, input_shape, output_shape);
%   
% Th and R are the radial co-ordinates of the input shape bin positions.
%
% Input_shape and output_shape define the transformation, as follows:
%
%   8   =   Square
% 	7   =   1:7 octagon
% 	6   =   2:6 octagon
% 	5   =   3:5 octagon
% 	4   =   4:4 (regular) octagon
%   1   =   Circle

a = find(th<0);                 % 
b(length(th)) = zeros;          % Need to change theta (0 to -pi) into (pi to 2*pi)
b(a) = 2*pi;                    %
th = th+b';                     %

shape_def = [shape_in, shape_out];
for i = 1:2    
    if shape_def(i) == 1        
        R{i} = ones([length(th), 1]);  %Rc(1:length(th),1) = 1;         % Radius of circle = 1, circumference = 2*pi.	    
    elseif shape_def(i) == 8       
        ss(1:length(th))        = 0;                    % 
        ss(th>(pi/4))     = 2;                    % ss is side indicator for square
        ss(th>((pi*3)/4)) = 4;                    %
        ss(th>((pi*5)/4)) = 6;                    %
        ss(th>((pi*7)/4)) = 0;                    %
        
        R{i} = (pi/4) ./ cos(((ss'.*pi)./4)-th);          % Make square perimeter
    else
		% Make Octagon, if necessary. (If input is not a circle). %
	
		ko = ((1+sqrt(2))*pi)./8;                       % radius of interior circle of regular octagon, for same circumference as above defined circle.
		
        n = shape_def(i);
		x = (n*pi)/16;                                  % Length of long side expressed in above terms (1 = radius of circle).
		
		regks = 1/(2+sqrt(8));                          % ks for regular octagon, for normalizing ks values to =1 for reg oct, see below.
        shortks = ((n/2) + sqrt(((8-n)^2)/2))*regks;    % ks is adjustment factor for irregular oct's - straight side is closer to centre than diagonal.
		longks  = (((8-n)/2) + sqrt((n^2)/2))*regks;    % These lines calculate absolute lengths (r's of interior circle(in morph box units)), but then express  
                                                        % that as a fraction of the value for regular oct, so shortks and longks for reg oct are =1.
		%%%%%%
		ths = atan(x/(2*shortks*ko));                   % Value of th for first corner.
		%%%%%%   Where are diagonal sides? (Where to use longks and shortks?) %%%%%%
		cd = [1,3,5,7].*(pi/4);                         % d is a row representing th values at the central points of the diagonals                      		
		stest = [cd(1)-th cd(2)-th cd(3)-th cd(4)-th];  % subtract each of d from th (looking for how close th to diag centres)		
		cl = min(abs(stest'));                          % For each value of th, how close to the closest diagonal	
		diag = find(cl<((pi/4)-ths));                   % Is closest th within the range defined by difference between first corner and first centre diagonal?	
		ks(1:length(th)) = shortks;                     % If so, ks = longks, ks = shortks otherwise.
		ks(diag) = longks;                              %   ...
		
        %%%%%%  Make side indicator, s. %%%%%%
		s(1:length(th)) = 0;
		s(th>ths)          = 1;
		s(th>(pi/2-ths))   = 2;
		s(th>(pi/2+ths))   = 3;
		s(th>(pi-ths))     = 4;
		s(th>(pi+ths))     = 5;
		s(th>(pi*3/2-ths)) = 6;
		s(th>(pi*3/2+ths)) = 7;
		s(th>(2*pi-ths))   = 0;
		%%%%%%                                                
		R{i} = (ko.*ks') ./ cos(((s'.*pi)./4)-th);        % Make octagon perimeter
    end 
end
% Transformation
rT = r .* (R{2}./R{1});      % Transformation: (radial distance from centre of bin) * ratio(trans_radius/input_radius), at that radial direction.




