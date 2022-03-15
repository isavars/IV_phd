function [KL]=dir_kldivergence(varValue,pVect1,pVect2,varargin)
% Kullback-Leibler or Jensen-Shannon divergence between two distributions.
%
%   KLDIV=dir_kldivergence(X,P1,P2) 
%
% returns the Kullback-Leibler divergence between two distributions 
% specified over the M variable values in vector X.  P1 is a length-M 
% vector of probabilities representing distribution 1, and P2 is a
% length-M vector of probabilities representing distribution 2.  Thus, the
% probability of value X(i) is P1(i) for distribution 1 and P2(i) for
% distribution 2.  The Kullback-Leibler divergence is given by:
%
%       KL(P1(x),P2(x)) = sum[P1(x).log(P1(x)/P2(x))]
%
% If X contains duplicate values, there will be an warning message, and these
% values will be treated as distinct values.  (I.e., the actual values do
% not enter into the computation, but the probabilities for the two
% duplicate values will be considered as probabilities corresponding to
% two unique values.)  The elements of probability vectors P1 and P2 must 
% each sum to 1 ± .00001.
%
% KLDIV(X,P1,P2,'sym') returns a symmetric variant of the Kullback-Leibler
% divergence, given by [KL(P1,P2)+KL(P2,P1)]/2.  See Johnson and Sinanovic
% (2001).
%
% KLDIV(X,P1,P2,'js') returns the Jensen-Shannon divergence, given by
% [KL(P1,Q)+KL(P2,Q)]/2, where Q = (P1+P2)/2.  See the Wikipedia article
% for "Kullback–Leibler divergence".  This is equal to 1/2 the so-called
% "Jeffrey divergence."  See Rubner et al. (2000).
%  
%   REFERENCES:
%   1) Cover, T.M. and J.A. Thomas. "Elements of Information Theory," Wiley, 
%      1991.
%   2) Johnson, D.H. and S. Sinanovic. "Symmetrizing the Kullback-Leibler 
%      distance." IEEE Transactions on Information Theory (Submitted).
%   3) Rubner, Y., Tomasi, C., and Guibas, L. J., 2000. "The Earth Mover's 
%      distance as a metric for image retrieval." International Journal of 
%      Computer Vision, 40(2): 99-121.
%   4) <a href="matlab:web('http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence','-browser')">Kullback–Leibler divergence</a>. Wikipedia, The Free Encyclopedia.
%
%   See also: MUTUALINFO, ENTROPY

if ~isequal(unique(varValue),sort(varValue)),
    warning('KLDIV:duplicates','X contains duplicate values. Treated as distinct values.')
end
if ~isequal(size(varValue),size(pVect1)) | ~isequal(size(varValue),size(pVect2)),
    error('All inputs must have same dimension.')
end
% Check probabilities sum to 1:
if (abs(sum(pVect1) - 1) > .00001) | (abs(sum(pVect2) - 1) > .00001),
    error('Probablities don''t sum to 1.')
end

if ~isempty(varargin),
    switch varargin{1},
        case 'js',
            logQvect = log2((pVect2+pVect1)/2);
            KL = .5 * (sum(pVect1.*(log2(pVect1)-logQvect)) + ...
                sum(pVect2.*(log2(pVect2)-logQvect)));

        case 'sym',
            KL1 = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));
            KL2 = sum(pVect2 .* (log2(pVect2)-log2(pVect1)));
            KL = (KL1+KL2)/2;
            
        otherwise
            error(['Last argument' ' "' varargin{1} '" ' 'not recognized.'])
    end
else
    zeroInd = pVect1==0 | pVect2==0;
    pVect1 = pVect1(~zeroInd);
    pVect2 = pVect2(~zeroInd);
    KL = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));
end


