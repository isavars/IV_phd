function [L_Ratio,IsoD,cellList] = spk_clusterisolation(amps,cellIndex)
% Calculate the quality of cluster separation on a tetrode, based on
% the L-Ratio of Schmitzer-Torbert & Redish, and the Isolation Distance
% of Harris (see Schmitzer-Torbert et al, 2005. Neuroscience 131).
%
%       [LR, IsoD, cellList] = spk_clusterisolation(amps,cellIndex)
%
% Inputs - amps is a [nSpike,4] array, containing the spike amplitude for
%          every spike on the tetrode, on each wire of the tetrode.
%        - cellIndex is an [nSpike,1] index vector to show which spikes belong to which cells.
%
% The results won't be meaningful unless all spikes are used, including
% 'Cell 0'. Cell 0 must be marked as 0 in the cellIndex vector.
%
% LR and IsoD will be vectors, one value for each cell (but not cell 0). cellList is a vector with the cell numbers.
%
% This function is an adaption of MClust3.4 (AD Redish) code.

nCells=length(unique(cellIndex)) - 1;

% Triode check %
active = find(sum(amps, 1) ~= 0); % Look for where all amps on wire are = 0
if length(active) == 3
    amps = amps(:,active);
elseif length(active) < 3
    disp('Won''t calculate cluster separation for stereotrodes or single wires');
    L_Ratio = NaN(nCells, 1);   IsoD = NaN(nCells, 1);
    return
end    

% Now get separation values, for each cluster at a time %
amps = double(amps);
L_Ratio = NaN(nCells, 1);
IsoD = NaN(nCells, 1);
cellList=setdiff(unique(cellIndex),0);
for ii=1:length(cellList)
    currCellInd = cellIndex==cellList(ii);
    noiseInd = ~currCellInd;
    if sum(currCellInd) <= size(amps,2) % N spikes must be greater than 
        continue                            % degrees of freedom for Mahalanobis
    end                                     % distance to work.
    
    % Mahalanobis distance %
    dists = mahal(amps,amps(currCellInd,:));
    
    % The following is taken from CLUSTERSEPARATION, part of MClust3.4.
    %------------------------------------
    % L_Ratio
    % take distances from mahal.m to be distributed as chi^2 with df
    % degrees of freedom [df = (# features)*(# channels)]
    % L_Extra is the sum of the probability that spikes from
    % outside the cluster (f_out) are in the probability density of the
    % cluster (f).
    df = size(amps,2);
    % Jadin suggested 1 - chi2pdf
    L_Extra = sum(1-chi2cdf(dists(noiseInd),df));
    L_Ratio(ii) = L_Extra/length(currCellInd);
    %-------------------------------------------
    % Ken Harris' cluster quality measure 
    % Another MClust3.4 function, included here as a sub-function.
    IsoD(ii) = ClusterQuality(amps,find(currCellInd),dists);    % ClusterQuality subfunction expects cell index as numerical vector
    %-------------------------------------------
    % End of MClust code
    
end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Isolation Distance Sub-function %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the function CLUSTERQUALITY, from the MClust3.4 toolbox.

function [HalfDist] = ClusterQuality(Fet, MySpikes, m)

% Measure of cluster quality
% [HalfDist] = ClusterQuality(Fet, MySpikes)
%
% see also FileQuality -  a wrapper function that runs this for every cluster
% in a given file.
%
% Inputs: Fet - N by D array of feature vectors (N spikes, D dimensional feature space)
% MySpikes: list of spikes corresponding to cell whose quality is to be evaluated.
% Res - Spike times
% BurstTimeWin - spikes within this time count as a burst
%
% make sure you only pass those features you want to use!
%
% Created by Ken Harris

% modified 18 Dec 02 ncst to accept m as an input, and a few formatting
% changes.

% find # of spikes in this cluster
if nargin < 3
	nSpikes = size(Fet,1);
else
	nSpikes = size(m,1);
end

nMySpikes = length(MySpikes);

% check there are enough spikes (but not more than half)
if length(MySpikes) < size(Fet,2) || length(MySpikes)>nSpikes/2
	HalfDist = 0;
	return
end

% mark other spikes
OtherSpikes = setdiff(1:nSpikes, MySpikes);

%%%%%%%%%%% compute mahalanobis distances %%%%%%%%%%%%%%%%%%%%%
if nargin < 3
	m = mahal(Fet, Fet(MySpikes,:));
end


mOther = m(OtherSpikes); % mahal dist of others

%fprintf('done mahalanobis calculation...');
% calculate point where mD of other spikes = n of this cell
if (nMySpikes < nSpikes/2)
	sorted = sort(mOther);
	HalfDist = sorted(nMySpikes);
else
	HalfDist = 0; % If there are more of this cell than every thing else, forget it.
end