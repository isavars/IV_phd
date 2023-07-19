function make_DS2_development_plots(electrodes)
%ds2 amp over age plot and other plots 

%TO DO - 1) dont include values from rats without a DS2 2) fix legends on 
% the per rat plots and gather mode data for elepos so they look nicer.  

    %load elePos 
    load(electrodes, 'elePos')
    %load useful parts from elePos
    age = elePos.age;
    rat_ID = elePos.rat_ID;
    %check that this is the best amplitude measure to use 
    ds2_amplitude = nanmax(abs(elePos.DS2_max_amplitude), [], 2);%make 1 max ds2_amplitude value per row of elePos
    ds2_rate = elePos.DS2_rate;

    %plots of DS2 amplitude over age 
    % Calculate mean and range per age
    [age_unique, ~, age_idx] = unique(age);
    mean_ds2_amplitude = splitapply(@nanmean, ds2_amplitude, age_idx);
    range_ds2_amplitude = splitapply(@(x) max(x) - min(x), ds2_amplitude, age_idx); %is this a normal standard error? 
 
    mean_ds2_rate = splitapply(@nanmean, ds2_rate, age_idx);
    range_ds2_rate = splitapply(@(x) max(x) - min(x), ds2_rate, age_idx); %is this a normal standard error? 


    % Plot mean with error bars
    figure;
    scatter(age_unique, mean_ds2_amplitude, 'filled', 'g');
    hold on;
    errorbar(age_unique, mean_ds2_amplitude, range_ds2_amplitude, 'LineStyle', 'none', 'Color', 'k');
    hold off;    
    % Adjust x-axis range
    xMin = min(age_unique) - 1;
    xMax = max(age_unique) + 1;
    xlim([xMin, xMax]);    
    % Add labels and title
    xlabel('Age');
    ylabel('Mean DS2 Max Amplitude');
    title('DS2 Max Amplitude over Age');

    %make another plot for rate over age 

    % Plot mean with error bars
    figure;
    scatter(age_unique, mean_ds2_rate, 'filled', 'r');
    hold on;
    errorbar(age_unique, mean_ds2_rate, range_ds2_rate, 'LineStyle', 'none', 'Color', 'k');
    hold off;    
    % Adjust x-axis range
    xMin = min(age_unique) - 1;
    xMax = max(age_unique) + 1;
    xlim([xMin, xMax]);    
    % Add labels and title
    xlabel('Age');
    ylabel('Mean DS2 firing rate');
    title('DS2 firing rate over Age');

    %figures per rat: 

    % Get unique rat IDs
    [unique_rat_IDs, ~, rat_idx] = unique(rat_ID);
    
    % Choose a colormap that will give a different color for each rat
    colors = jet(length(unique_rat_IDs));
    
    % Create a new figure for the amplitude over age plot
    figure;
    
    % A cell array to hold the labels for the legend
    labels = cell(length(unique_rat_IDs), 1);
    
    % Loop over the rat IDs
    for i = 1:length(unique_rat_IDs)
        % Select the data for this rat
        this_rat_idx = rat_idx == i;
        this_rat_age = age(this_rat_idx);
        this_rat_amplitude = ds2_amplitude(this_rat_idx);
    
        % For age 40, replace all values with their mean
        is_age_40 = this_rat_age == 40;
        if any(is_age_40)
            this_rat_amplitude(is_age_40) = mean(this_rat_amplitude(is_age_40));
        end
    
        % Sort the data by age
        [this_rat_age, idx] = sort(this_rat_age);
        this_rat_amplitude = this_rat_amplitude(idx);
    
        % Create a scatter plot for this rat's data
        scatter(this_rat_age, this_rat_amplitude, 'filled', 'MarkerFaceColor', colors(i, :));
        hold on;
    
        % Add a line connecting the dots
        plot(this_rat_age, this_rat_amplitude, 'Color', colors(i, :));
    
        % Add this rat's ID to the labels for the legend
        labels{i} = unique_rat_IDs{i};
    end
    
    % Create the legend
    legend(labels);
    
    % Add labels and title
    xlabel('Age');
    ylabel('DS2 Max Amplitude');
    title('DS2 Max Amplitude over Age per Rat');
    hold off;
    
    % Create a new figure for the rate over age plot
    figure;
    
    % A cell array to hold the labels for the legend
    labels = cell(length(unique_rat_IDs), 1);
    
    % Loop over the rat IDs
    for i = 1:length(unique_rat_IDs)
        % Select the data for this rat
        this_rat_idx = rat_idx == i;
        this_rat_age = age(this_rat_idx);
        this_rat_rate = ds2_rate(this_rat_idx);
    
        % For age 40, replace all values with their mean
        is_age_40 = this_rat_age == 40;
        if any(is_age_40)
            this_rat_rate(is_age_40) = mean(this_rat_rate(is_age_40));
        end
    
        % Sort the data by age
        [this_rat_age, idx] = sort(this_rat_age);
        this_rat_rate = this_rat_rate(idx);
    
        % Create a scatter plot for this rat's data
        scatter(this_rat_age, this_rat_rate, 'filled', 'MarkerFaceColor', colors(i, :));
        hold on;
    
        % Add a line connecting the dots
        plot(this_rat_age, this_rat_rate, 'Color', colors(i, :));
    
        % Add this rat's ID to the labels for the legend
        labels{i} = unique_rat_IDs{i};
    end
    
    % Create the legend
    legend(labels);
    
    % Add labels and title
    xlabel('Age');
    ylabel('DS2 firing rate');
    title('DS2 firing rate over Age per Rat');
    hold off;



end