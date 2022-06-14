function rows = selectrows(T, key)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
keyFields = fieldnames(key);
varNames = T.Properties.VariableNames;
selIdx = true(height(T), 1);
for ii = 1:length(varNames)
    if ismember(varNames{ii}, keyFields)
        column = T.(varNames{ii});
        value = key.(varNames{ii});
        if isnumeric(column)
            selIdx = selIdx & (column == value);
        else % should be a cell array of strings
            selIdx = selIdx & arrayfun(@(x) strcmp(x, value), column);
        end
    end
end
rows = T(selIdx,:);
end

