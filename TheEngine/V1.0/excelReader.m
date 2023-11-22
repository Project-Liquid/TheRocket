function excelValue = excelReader(name, excelTable)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    rowIndex = find(strcmp(excelTable{:,1}, name));
    rowValue = excelTable(rowIndex, 2);
    excelValue = rowValue{1,1};
end

