rows_to_clear = [1,2,14, 15, 26, 27, 28, 29, 30, 31, 32, 34,35,36,37]; % requires visual inspection of DS2 figures for each trial. 

for i = rows_to_clear
    elePos.DS2_spike_times{i} = {};
    elePos.DS2_labels(i,:) = nan;
    elePos.DS2_max_amplitude(i,:) = nan;
    elePos.DS2_mean_amplitude(i,:) = nan;
    elePos.DS2_peak_to_trough_amplitude(i,:) = nan;
    elePos.DS2_rate(i,:) = nan;
    elePos.DS2_slope(i,:) = nan;
end

    
  