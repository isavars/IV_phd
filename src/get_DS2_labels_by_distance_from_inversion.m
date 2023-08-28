function DS2_orientations = get_DS2_labels_by_distance_from_inversion(data, electrodes)
%needs to create proper ds2 inverison metric with distance from inversion
%will need to read in channel depths from maps in office computer - add
%these to elePos 

%for elePos or a different metadata file to be used. 
load (electrodes, 'elePos');
%load spatial Data 
load (data, 'spatData');


end