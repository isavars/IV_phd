function get_DS2_cell_cofiring(data, electrodes)
% need to think about this a bit more mathematically but also check S&B
% methods - basically this function needs to read in the spike times of
% cells and compare them to the times the ds2 fire - maybe this produces
% histograms like the ones in S&B - 
    % Peri-event time histogram of GC (left) and MC (right) spikes aligned to 
    % the time of the peak of DS2 (mean ± s.d.). 
    % Peak firing ratio was calculated by ..
    % I think what this means is i have to take 1 second time windows
    % around each DS2 and plot the ratio of the maximum firing rate of over 
    % the baseline firing rate of a cell for each time unit in that window 
    % - then plot all the spike windows overlayed for that group (gc or mc)

    % the relevant result is that GCs are more coincident with ds2 than
    % mossy but they both are and for CA3 i guess it should go down (if you
    % can get CA3 on a rat that also has DS2 for sure. 

    %the other thing for classifing ca3 could just be is there or isnt
    %there a DS2 on the tetrode the cell is recorded on? like if they were
    %detected for a bunch of days and then you turned again but have no ds2
    %thats ca3 probably - hard to say per individual terode though 

% this function should be used in two ways: 
%   - pre classification to be used in classifying CA3 cells out 
%   - post classification to provide results about the DS2 activity vs 
%     the firing of cells different kinds - in development - appropriate
%     stats need to be conducted to compare pre and post wean 
%   Inputs - needs data from spike data for cell identities and elepos for
%   ds2 spike times. Outputs - will output for every row of spatData (so
%   per cell) if the cells firing co-incided with the ds2 

end