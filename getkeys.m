function keys = getkeys(conditions, selections, varargin)
%GETKEYS Summary of this function goes here
%   Detailed explanation goes here
P = inputParser;
addRequired(P, 'conditions', @istable)
addRequired(P, 'selections', @isstruct)
addOptional(P, 'forces', table(), @istable)
addOptional(P, 'psths', table(), @istable)
addOptional(P, 'spikes', table(), @istable)
parse(P, conditions, selections, varargin{:})

% join input data with optional tables
data = conditions;

baseJoinKeys = {'experiment', 'session_index', 'condition_block', 'condition_index'};

if ~isempty(P.Results.forces)
    data = innerjoin(data, P.Results.forces, 'LeftKeys', baseJoinKeys, 'RightKeys', baseJoinKeys); 
end

if ~isempty(P.Results.spikes)
    if ismember('trial_number', data.Properties.VariableNames)
        joinKeys = [baseJoinKeys, {'trial_number'}];
    else
        joinKeys = baseJoinKeys;
    end
    data = innerjoin(data, P.Results.spikes, 'LeftKeys', joinKeys, 'RightKeys', joinKeys);
end

if ~isempty(P.Results.psths)
    if ismember('motor_unit_index', data.Properties.VariableNames)
        joinKeys = [baseJoinKeys, {'motor_unit_index'}];
    else
        joinKeys = baseJoinKeys;
    end
    data = innerjoin(data, P.Results.psths, 'LeftKeys', joinKeys, 'RightKeys', joinKeys);
end

% identify scalar data fields
varNames = data.Properties.VariableNames;
isScalar = false(1, length(varNames));
for ii = 1:length(varNames)
    value = table2cell(data(1, varNames{ii}));
    isScalar(ii) = (isscalar(value{1}) || ischar(value{1}));
end
scalarFields = varNames(isScalar);

% update selection with all valid scalar fields
selections = selections(:);
keys = [];
for ii = 1:length(selections)
    dataSel = selectrows(data, selections(ii));
    if height(dataSel) ~= 1
        error('Insufficient selectors for input arguments.')
    end
    
    key = struct();
    for jj = 1:length(scalarFields)
        value = table2cell(dataSel(1, scalarFields{jj}));
        key.(scalarFields{jj}) = value{1};
    end
    keys = [keys; key];
end

end

