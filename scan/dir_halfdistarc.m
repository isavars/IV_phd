function [rtn]=dir_halfdistarc(map)
% Calculates the size of the arc around the circular mean that contains
% half of the distribution of firing rate.
% See Solstad et al, 2008 (Science).

if map==0
    rtn=nan;   return
end

% Circular Mean %
binSize=360/length(map);
binAngles=(0:binSize:359) + binSize/2;
circMean= circ_rad2ang(circ_mean( circ_ang2rad(binAngles)', map ));
if circMean<0;   circMean=circMean+360;   end

% Realign dir map so that circ mean bin is map(1)
meanBin=ceil((circMean+1)./binSize);
shiftMap=map([meanBin:length(map) 1:meanBin-1]);

% Step through increasing arcs until half Dist threshold is crossed %
halfDistr=sum(map)/2;
indCount=0;
arc=0;
while sum(arc)<halfDistr
    indCount=indCount+1;
    arc=[shiftMap(1:indCount); shiftMap(end-indCount+2:end)];
end
rtn=length(arc)*binSize;

