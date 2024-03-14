function params = excelReader(excelPath, variablesList)
    % ExcelPull reads data from an Excel file and creates variables dynamically
    
    % Read the Excel file into a table
    data = readtable(excelPath, 'PreserveVariableNames', true);
    numCols = width(data);
    params = struct();

    % Iterate through the variables list
    for i = 1:length(variablesList)
        variableName = variablesList{i};
        found = false;

        % Iterate over all columns (except the last one)
        for col = 1:numCols-1
            % Search in the current column
            rowIndex = find(strcmp(data{:, col}, variableName));

            if ~isempty(rowIndex)
                % Value is assumed to be in the next column
                value = data{rowIndex, col+1};
                params.(variableName) = value;
                found = true;
                break; % Exit the column loop as the variable is found
            end
        end

        if ~found
            disp(['Variable ' variableName ' not found in Excel file.']);
        end
    end
end
