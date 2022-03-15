function R = deg2rad(D)

if nargin==0
    error('Incorrect number of arguments')
elseif ~isreal(D)
    warning('Imaginary parts of complex ANGLE argument ignored')
    D = real(D);
end

% D = R*180/pi;
R = D*pi/180;

% ----------------------------------------------------------------------------------------------------------------------
