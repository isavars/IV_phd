function D = rad2deg(R)

if nargin==0
    error('Incorrect number of arguments')
elseif ~isreal(R)
    warning('Imaginary parts of complex ANGLE argument ignored')
    R = real(R);
end

D = R*180/pi;

% ----------------------------------------------------------------------------------------------------------------------
