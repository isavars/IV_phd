function [edgeLengths, edgeRegularity, nnLengths, nnRegularity]=map_fieldspacing(com)
% Scores of evenness of field spacing.
% Input: array of centre-of-mass field centres.

if size(com,1)<3
    edgeLengths=nan; edgeRegularity=nan; nnLengths=nan; nnRegularity=nan;
    return
end
    
%%%%%%%%%%%%% Part 1: Nearest-neighbour analyses %%%%%%%%%%%%%%%%%%%
interFieldDist=nan(size(com,1));
for ii=1:size(com,1)
    for jj=1:size(com,1)
        if ii~=jj
            interFieldDist(ii,jj)=sqrt( (com(ii,1)-com(jj,1))^2 + (com(ii,2)-com(jj,2))^2 );
        end
    end
end
nnLengths=nanmin(interFieldDist);
nnRegularity=std(nnLengths)/mean(nnLengths);

%%%%%%%%%%%%% Part 2: Delauney triangulation analyses %%%%%%%%%%%%%%
%%%%%%%% Copied from CB DELAUNAY_GRID_ANALYSIS, 04/02/10 %%%%%%%%%%%
try
    tri=delaunay(com(:,1), com(:,2), {'Qt','Qbb','Qc','Qz'});
catch
    edgeLengths=nan; edgeRegularity=nan;
    return
end
% Remove sides on the convex hull %
% K = convhull(com(:,1), com(:,2));
% for ii=1:length(com)
    

% --- GENERATE MEASURES BASED ON THE TRIANGLES WE'VE FOUND ------------------------------
%First need to find unique edges - triangles reuse edges so dont want to count them twice
shiftOne=[2,3,1];
firstVert=[];
secondVert=[];
for nn=1:size(tri,1)
    for kk=1:3
        currVert=sort([tri(nn,kk), tri(nn,shiftOne(kk))]);
        if isempty(firstVert)
            firstVert=currVert(1);
            secondVert=currVert(2);
        elseif ~any((firstVert==currVert(1)) .* (secondVert==currVert(2)))
            firstVert(end+1)=currVert(1);
            secondVert(end+1)=currVert(2);
        end
    end
end
firstCoord=com(firstVert,:); %x,y pair with origin top left
secondCoord=com(secondVert,:);
%Get distance between peaks and angle between them - note I've flipped the y axis to take
%account of the fact that origin is top left
[angle, edgeLengths]=cart2pol(firstCoord(:,1)-secondCoord(:,1), secondCoord(:,2)-firstCoord(:,2));
edgeLengths=abs(edgeLengths);
edgeRegularity=std(edgeLengths)/mean(edgeLengths);
angle=mod(rad2deg(angle),360);
angleMod60=mod(angle,60);
meanAngle=mod(rad2deg(circ_mean(deg2rad(angleMod60.*6)))/6, 60);


