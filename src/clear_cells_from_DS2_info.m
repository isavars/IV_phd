rows_to_clear = [1,2]

for i = rows_to_clear
    elePos.DS2_spike_times{i} ={}
    elePos.DS2_labels(i,:) = nan
end

    DS2_info.DS2_max_amplitude(1:32) = nan;
    DS2_info.DS2_mean_amplitude(1:32) = nan;
    DS2_info.DS2_rate(1:32) = nan;
    DS2_info.DS2_peak_to_trough_amplitude(1:32) = nan;
    DS2_info.DS2_slope(1:32) = nan;    
  