% Function for jittering data points
function jitteredData = jitter(data)
    jitterRange = 0.5;  % Adjust as needed
    jitteredData = data + (rand(size(data))-0.5) * jitterRange;
end