function [p, d, converged] = pxd_mtint(spikes, times, max_iter, accuracy, tol)
% [p, d, converged] = mtint_pxd(spikes, times, max_iter, accuracy, tol);
% plots the fields and polar_plot resulting from a maximum likelihood factorial
% model of influences of place (p) and direction (d) on cell firing, assuming Poisson noise.
% Standard (spikes/dwell_time) plots are also shown for comparison.
% See Burgess, Cacucci, Lever, O'Keefe Hippocampus (2005) for details of procedure, although
% slightly different binning and smoothing of data is used in this reference.
%
% 'spikes' and 'times' are matrices with dimensions: d_bins y_bins x_bins.
%
% The algorithm is run until the fractional change in the likelihood of the data
% is less than 'accuracy' (e.g. 0.00001, i.e. until 'convergence') or for max_iter iterations
% (e.g. 30, i.e. 'non-convergence').
%
% The value tol (e.g. 0.1) replaces expected firing rate values less than tol (usually
% due to unvisited states) to avoid log(0) divergence in loglikelhood, and to avoid division
% by zero in p_estimate and d_estimate.
%
% Output: p is a matrix with dimensions: y_bins x_bins (y increases downwards),
% d is a vector with dimension d_bins.
%
% Use e.g.:
% [p, d, converged] = pxd_mtint(spikes, times, 30, 0.000001, 0.1);

converged = 0;
if ( sum( size(spikes) == size(times) ) ~= 3 )
    fprintf(1, ' matrices spikes and times have different dimensions - spikes: %d %d %d times: %d %d %d\n',...
        size(spikes), size(times));
else
    [nd ny nx] = size( spikes);
    fit = 1;
%     fprintf(1, ' loglikelihood:');
    for iter = 1:max_iter
        if( iter==1 )
            % Initial guess
            p = ones( ny, nx);
            d = ones( nd );
        else
            p = p_estimate(d, spikes, times, tol);
        end
        d = d_estimate(p, spikes, times, tol);
        prev_fit = fit;
        %
        % NB log likelihhod is a negative (the smaller in magnitude the better)
        %
        fit = loglikelihood(p, d, spikes, times, tol);
        if( abs(prev_fit - fit) < -accuracy*fit )
%             fprintf(1, ' converged, loglikelihood: %f\n', fit);
            converged = 1;
            break;
        end
%         fprintf(1, ' %f', fit);
    end

    if( converged == 0 )
        fprintf(1, ' Did not converge, %d spikes, nd=%d ny=%d nx=%d.\n', sum(sum(sum(spikes))), size(spikes));
    else
        % for direct comparison with regular fields, normalise area under p and d separately to
        % equal the tot no spikes / tot dwell time.
        tot_spikes = sum(sum(sum(spikes)));
        pred_spikes = sum(sum(p.*(squeeze(sum(times,1)))));
        p = p.*(tot_spikes/pred_spikes);
        pred_spikes = sum(d.*(squeeze(sum(sum(times,2),3))));
        d = d.*(tot_spikes/pred_spikes);
    end
end