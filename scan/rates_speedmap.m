function [rtn] = rates_speedmap(y,x,dir,speed,window,bins)
% Makes a 'rate map' showing the mean speed in each spatial bin.
%
%   [rtn] = rates_speedmap(y,x,speed,window,bins)
%
% 'x', 'y' 'dir' and 'speed' are vectors of these values for each pos sample, from the trial data. x,y in pix, dir in degrees and speed in cm/s.
% 'window' sets the max extent of a 2D place map, format [yExtent xExtent], in pix. Lower limits for map are always [0 0].
% 'bins' is in the format [yBin xBin dirBin]. Currently, can only do EITHER place OR dir (not PxD), so either dirBin should be set to 360 (for place)
% or to <360 for a dir map. If dirBin<360, the x and y bins are ignored.


if bins(3)==360
    binnedY = ceil( double(y) ./bins(1) );
    binnedX = ceil( double(x) ./bins(2) );
    binnedLinInd = sub2ind(ceil(window./bins(1:2)), binnedY, binnedX);
else
    binnedLinInd = ceil( double(dir) ./bins(3) );
end

binList = unique(binnedLinInd);
dwellHist = histc(binnedLinInd,0.5:1:(max(binnedLinInd)+0.5));

sampSortByBin = nan(max(dwellHist),prod(   ceil(window./bins(1:2))   ) );
for ii=1:length(binList)
    binIndTemp = binnedLinInd==binList(ii);
    sampSortByBin(   1:sum(binIndTemp),   binList(ii)  ) = speed(binIndTemp);
end
meanSpeedLinearMap = nanmean(sampSortByBin);

if bins(3)==360
    rtn = reshape(meanSpeedLinearMap,ceil(window./bins(1:2)));
else
    rtn = meanSpeedLinearMap;
end

    
%     binList = unique(binnedLinInd);
%     errorMap = nan(25,25);
%     posMap = nan(25,25);
%     for jj=1:length(binList)
%         binIndTemp = binnedLinInd==binList(jj);
%         if any(     strcmp(errorMetricType,{'correctPredDistribution','none (plot pos)'})     )
%             errorMap(   binList(jj)   ) =  sum(     errorMetric( binIndTemp )    );  
%         else
%             errorMap(   binList(jj)   ) = nanmean(     errorMetric( binIndTemp )    );
%         end
%         posMap(   binList(jj)    ) = sum(  binIndTemp  );
%     end