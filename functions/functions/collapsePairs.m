function collapsedStruct = collapsePairs(data, pairs)
% COLLAPSEPAIRS Collapse experiments, types, and ori tables into a single struct
%
% Inputs:
%   data  - struct with fields:
%           .exp  : cell array of experiment names
%           .type : cell array of types
%           .ori  : cell array of tables/matrices
%   pairs - cell array of vectors specifying which experiments to collapse
%
% Output:
%   collapsedStruct - 1x1 struct with fields:
%                     .exp, .type, .ori (each is a cell array of length numel(pairs))

nPairs = numel(pairs);

expField  = cell(1,nPairs);
typeField = cell(1,nPairs);
oriField  = cell(1,nPairs);

for ii = 1:nPairs
    idx = pairs{ii};

    % Merge experiment names
    expField{ii}  = strjoin(data.exp(idx), '_');

    % Merge types
    typeField{ii} = strjoin(data.type(idx), '_');

    % Merge ori tables/matrices vertically
    oriField{ii}  = vertcat(data.ori{idx});
end

% Build the single struct
collapsedStruct.exp  = expField;
collapsedStruct.type = typeField;
collapsedStruct.ori  = oriField;

end
