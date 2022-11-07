function [tetrode_new,num_spikes] =  cut_trial_to_time(tetrode_data,trial_duration,num_spikes)
%% Cut trial to correct trial length (usually 10 mins)
% INPUT:
% time_mat: array of timestamps (output from concat_trial func)
% complete_trial: nx16 matrix of voltages for trial (output from concat_trial func)
% trial_length: trial length in mins
% 
% OUTPUT:
% time_new: array of timestamps to trial end
% complete_trial_new: nx16 matrix of voltages to trial end

% trial_duration = trialInfo.trial_duration;
% trial_duration = trial_duration * 60; % seconds
length_microseconds  = trial_duration * 1e6; % microseconds

time_excess = find(tetrode_data(:,1) > length_microseconds,1); % time in excess of duration

tetrode_new = tetrode_data;
tetrode_new(time_excess:end,:) = [];

num_spikes = num_spikes - numel(tetrode_data(time_excess:end,1))/4;

end

