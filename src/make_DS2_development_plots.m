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

    
    % Linear regression for mean_ds2_rate against age_unique
    lm = fitlm(age_unique, mean_ds2_rate, 'linear');
    
    % Display the results of the simple regression
    disp(lm);
    
    % Create a scatter plot for mean_ds2_rate vs age_unique with the regression line
    figure;
    scatter(age_unique, mean_ds2_rate, 'filled', 'r');
    hold on;
    plot(age_unique, lm.Fitted, '-b', 'LineWidth', 2);
    xlabel('Age');
    ylabel('Mean DS2 firing rate');
    title('Regression of DS2 firing rate over Age (aggregated)');
    hold off;


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
    
%     % Create the legend
%     legend(labels);
    
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
    %legend(labels);
    
    % Add labels and title
    xlabel('Age');
    ylabel('DS2 firing rate');
    title('DS2 firing rate over Age per Rat');
    hold off;
    
    %stats 

    % Create a table with your data
    data_table = table(age, rat_ID, ds2_amplitude, 'VariableNames', {'Age', 'Rat', 'Amplitude'});
    
    % Fit the linear mixed-effects model
    lme = fitlme(data_table, 'Amplitude ~ 1 + Age + (1 + Age | Rat)'); %this includes random intercepts which means the baseline amplitudes of each rat are considered to be different. 
    
    % Display the results
    disp(lme);

    % Create a table with your data for rate 
    data_table = table(age, rat_ID, ds2_rate, 'VariableNames', {'Age', 'Rat', 'Rate'});
    
    % Fit the linear mixed-effects model
    lme = fitlme(data_table, 'Rate ~ 1 + Age + (1 + Age | Rat)'); %this includes random intercepts which means the baseline amplitudes of each rat are considered to be different. 
    
    % Display the results
    disp(lme);

    %checking to see data is normally distributed

    % Obtain residuals
    residuals = lme.Residuals.Raw;
    
    % Create a histogram to visualize the distribution of the residuals
    figure;
    histogram(residuals, 'Normalization', 'pdf');
    hold on;
    
    % Compare the histogram with a normal distribution
    x = linspace(min(residuals), max(residuals), 100);
    y = normpdf(x, mean(residuals), std(residuals));
    plot(x, y, 'LineWidth', 2);
    
    % Add labels and title
    xlabel('Residuals');
    ylabel('Probability Density');
    title('Distribution of Residuals');
    hold off;
    
    % Create a scatter plot of residuals vs. fitted values
    figure;
    scatter(lme.fitted, residuals);
    
    % Add labels and title
    xlabel('Fitted values');
    ylabel('Residuals');
    title('Residuals vs. Fitted Values');

    predicted_values = predict(lme);
    % 1. Plot residuals against predicted values (for homoscedasticity)
    figure;
    scatter(predicted_values, residuals);
    xlabel('Predicted Values');
    ylabel('Residuals');
    title('Residuals vs Predicted Values');
    
    % If the plot shows a funnel shape (residuals fanning out as the predicted 
    % values increase), it's a sign that the assumption of homoscedasticity might be violated.
    
    % 2. Q-Q plot of residuals (for normality)
    figure;
    qqplot(residuals);
    title('Q-Q Plot of Residuals');
    
    % If the points in this plot lie along the reference line, it's a good sign 
    % that the residuals are normally distributed. If they deviate substantially, 
    % the assumption of normality might be violated.
    
 
    % If the autocorrelations are near zero for all lags, it's a good sign that 
    % the assumption of independence is satisfied. If there are one or more lags 
    % with autocorrelations that are substantially different from zero, the 
    % assumption of independence might be violated.

end

