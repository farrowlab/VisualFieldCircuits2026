function Coordinates = load_cell_coordinates(dataFolder, expName, areasToRemove)
%LOAD_CELL_COORDINATES Load the labelled cells of one tracing experiment.
%   Coordinates = LOAD_CELL_COORDINATES(dataFolder, expName, areasToRemove)
%   loads the table 'Coordinates' (one row per labelled cell) from
%   <dataFolder>/<expName>.mat, makes the AP coordinates positive and
%   removes cells in areasToRemove (see remove_areas).

S = load(fullfile(dataFolder, [expName '.mat']), 'Coordinates');
Coordinates = S.Coordinates;

% Make AP coordinates positive in all brains
if median(Coordinates.AP_location) <= 0
    Coordinates.AP_location = -Coordinates.AP_location;
end

Coordinates = remove_areas(Coordinates, areasToRemove);
end
