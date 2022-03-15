function [dotProdMean,dotProd]=map_popvect(mapsA,mapsB)
% Population vector comparison of ensemble firing in two different trials
%
%       [C,D]=map_popvect(mapsA,mapsB)
%
% Where C is the mean dot product between ensemble firing rate vectors from each
% place bin, between trials A and B. The dot product is normalised to vector
% magnitude, so C represents the mean cosine between vectors, ranging from 
% 1 (max similarity) to 0. Note that this measure therefore ignores absolute firing rates.
%
% D is an array of indivdual dot products from each bin, with the same spatial layout as the rate maps.
%
% mapsA and mapsB are cell arrays of the rate maps from each trial, in the format {1,nCell}. Each 
% map within the cell array is an (nBin,nBin) 2D matrix. Non-visited bins should be set to NaN.

stackA = cat(3,mapsA{:});
stackB = cat(3,mapsB{:});

visMaskA = sum(isnan(stackA),3) == 0;
visMaskB = sum(isnan(stackB),3) == 0;
visMask = visMaskA & visMaskB;

dotProd=nan(size(stackA,1),size(stackA,2));
for ii=1:size(dotProd,1)
    for jj=1:size(dotProd,2)
        if visMask(ii,jj)
            a = double(squeeze( stackA(ii,jj,:) ));
            b = double(squeeze( stackB(ii,jj,:) ));
            dotProd(ii,jj) = dot(a,b) / (norm(a)*norm(b));
        end
    end
end
dotProdMean=nanmean(dotProd(:));