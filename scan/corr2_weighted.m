function [r_out] = corr2_weighted( X, Y, W )
% Weighted correlation between two vectors or matrices of the same size.
%
%       [r_out] = corr2_weighted( X, Y, W )
%
% Like matlab CORR2, but can include a weight vector/matrix, W, to calculate
% the weighted version of the Pearson'r correlation between X and Y.
% It can also handle nans in the input - any points which are nan for X, Y, or W
% are excluded from the correlation.

isValid = ~isnan(X) & ~isnan(Y) & ~isnan(W);
X = X(isValid);
Y = Y(isValid);
W = W(isValid);

W = W ./ sum(W);   % Making W sum to 1 means we don't have to write the denominator in every eq.

weightedMeanX = sum( X .* W );
weightedMeanY = sum( Y .* W );

weightedVarX  = sum( W .* ((X - weightedMeanX).^2));
weightedVarY  = sum( W .* ((Y - weightedMeanY).^2));

weightedCov   = sum( W .* (X - weightedMeanX) .* (Y - weightedMeanY) );

r_out         = weightedCov / sqrt( weightedVarX * weightedVarY );
