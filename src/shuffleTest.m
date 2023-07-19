function spatData = shuffleTest(spatData, numShuffles)

    % iterate over rows/neurons in spatData
    for i = 1:height(spatData)
        % get the firing rate map, and original spatial information score for the current neuron
        originalFRMap = spatData.rMap{i};
        originalScore = spatData.SI_spat(i);

        % compute the total time length and the number of bins for shuffling
        [num_rows, num_cols] = size(originalFRMap);
        num_bins = num_rows * num_cols;
        reshapedFRMap = reshape(originalFRMap, [], num_bins);

        % initialize the array to store the scores from shuffled data
        shuffledScores = zeros(numShuffles, 1);

        % run shuffling procedure for numShuffles iterations
        for j = 1:numShuffles
            % shuffle the order of bins in the firing rate map
            shuffledFRMap = reshapedFRMap(:, randperm(num_bins));
            % flatten the shuffled map back into the original shape
            shuffledFRMap = reshape(shuffledFRMap, size(originalFRMap));

            % generate the times array for map_skaggsinfo
            times = ones(size(shuffledFRMap));

            % compute spatial information score with shuffled firing rate map
            [~, bits_per_sec] = map_skaggsinfo(shuffledFRMap, times);
            shuffledScores(j) = bits_per_sec;
        end

        % test if all shuffled scores are less than the original score
        if all(shuffledScores < originalScore)
            spatData.sig_SI(i) = true;
        else
            spatData.sig_SI(i) = false;
        end
    end
end

% Your function to calculate spatial information score
function [bits_per_spike,bits_per_sec] = map_skaggsinfo(rates, times)
    if max(max(rates)) <= 0
        bits_per_spike = NaN;
        bits_per_sec = NaN;
        return
    end

    if (size( rates ) ~= size(times))
        error('Size mismatch between spikes and position inputs to skaggs_info');
    end

    rates(isnan(rates)) = 0;   % remove background
    times(isnan(times)) = 0;

    if(size( rates ) > 1)
        % turn arrays into column vectors
        rates = reshape( rates, prod(size(rates)), 1);
        times = reshape( times, prod(size(times)), 1);
    end

    duration = sum(times);
    mean_rate = sum(rates.*times)./duration;

    p_x = times./duration;
    p_r = rates./mean_rate;                   
    dum = p_x.*rates;        
    ind = find( dum > 0 );   
    bits_per_sec = sum(dum(ind).*log2(p_r(ind)));   
    bits_per_spike = bits_per_sec/mean_rate;
end
