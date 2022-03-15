function [meanR, meanDir]=dir_rayleighvector(map)
% Rayleigh vector (mean resultant length) derived from directional firing 
% rate map, i.e. vector sum of map weighted by directional bin rate value.
%       [rtn]=dir_rayleighvector(map)

if isempty(map);   meanR=nan; meanDir=nan;   return;   end
if sum(map)==0;   meanR=nan; meanDir=nan;   return;   end

binSize=(pi*2)/length(map);
binAngles=(0:binSize:( (359.5/360)*2*pi ))' + binSize/2;
binWeights=map./(nanmax(map));
%calculate rayleigh vector
ral_vect = nansum(binWeights .* exp(1i*binAngles)); %complex 
%grab mean length and mean angle of vector
meanR = abs(ral_vect)/nansum(binWeights); %real part normalised by weights 
meanDir = angle(ral_vect); %complex part is angle
% S=nansum( sin(binAngles).*binWeights );
% C=nansum( cos(binAngles).*binWeights );
% 
% R=sqrt(S^2+C^2);
% % phi=((atan2(S,C)+pi)/pi*2)*360;
% meanR=R/nansum(binWeights);


end