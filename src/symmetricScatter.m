function symmetricHistogramScatter(data, binWidth, midpoint)
    if nargin < 3
        midpoint = 0;
    end

    if nargin < 2
        binWidth = 1; % Default bin width
    end
    
    % Calculate the histogram
    [counts, edges] = histcounts(data, 'BinWidth', binWidth);
    binCenters = edges(1:end-1) + diff(edges) / 2;

    figure;
    hold on;

    for i = 1:length(counts)
        numPoints = counts(i);
        
        % If we have an even number, we spread symmetrically with an equal number of points on either side of the midpoint.
        % If odd, the midpoint will have one point.
        if mod(numPoints, 2) == 0
            half = numPoints / 2;
            xvals = [-half + 0.5 : 1 : half - 0.5] + midpoint;
        else
            half = (numPoints - 1) / 2;
            xvals = [-half:1:half] + midpoint;
        end
        
        % Y-values are the bin centers
        yvals = repmat(binCenters(i), 1, numPoints);
        
        scatter(xvals, yvals, 'filled');
    end
    
    % Aesthetics
    xlabel('Symmetric Spread around Midpoint');
    ylabel('Data Value');
    title('Symmetric Histogram Scatter');
    grid on;

    hold off;
end
