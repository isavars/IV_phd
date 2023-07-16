%adaptable inputs
old_string = '190525_P20';
new_string = '190525_P22';
table = spatData;

% Get the number of rows in the table
numRows = size(table, 1);

% Iterate over each row
for row = 1:numRows
    % Check if cellID contains old_string 
    if contains(table.cellID{row}, old_string)
        % Replace old_string with new_string 
        table.cellID{row} = strrep(table.cellID{row}, old_string, new_string);
    end
    
    % Check if dataset contains old_string
    if contains(table.dataset{row}, old_string)
        % Replace old_string with new_string
        table.dataset{row} = strrep(table.dataset{row}, old_string, new_string);
    end
end

spatData = table; 