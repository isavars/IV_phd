function [com labelRm, nFields]=map_findfields(rm,fieldCutOff,minSize)
% Find the fields centres-of-mass 
% Taken from CB delaunay_grid_analysis, 04/02/10
labelRm=rm>fieldCutOff;
labelRm=bwlabel(labelRm,4);
stats=regionprops(labelRm, rm, 'area', 'weightedcentroid' );
ind = vertcat(stats.Area)<minSize;
invalidAreas=find(ind);
stats=stats(~ind);
com=vertcat(stats.WeightedCentroid); %Centre of mass of each peaks weighted by firing rate
% Remove too small fields and re-label %
for ii=1:length(invalidAreas)
    labelRm(labelRm==invalidAreas(ii))=0;
end
labelRm=bwlabel( labelRm>0 );
nFields=length(unique( labelRm )) - 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 09/10/10: Function below has been decommisioned as an independent
% function. Needs to be integrated into this function.

% function [com labelRm]=map_findpeaks(rm,peakThr,closePeakFilter)
% % Find peaks ie local maxima above a threshold. Can also specify a minimum separation for peaks.
% rm(isnan(rm))=0;
% peaks = imregionalmax(rm);
% peaks(rm<peakThr)=0;
% peaks = imdilate(peaks, closePeakFilter); 
% peaks = imfill(peaks,'holes');
% labelRm=bwlabel(peaks,4);
% stats=regionprops(labelRm, 'Centroid');
% com=vertcat(stats.Centroid);