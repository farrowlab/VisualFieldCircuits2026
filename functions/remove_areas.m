function Coordinates = remove_areas(Coordinates, areasToRemove)
%REMOVE_AREAS Remove labelled cells located in the given areas.
%   Coordinates = REMOVE_AREAS(Coordinates, areasToRemove) removes rows
%   whose area acronym (column 2) starts with any of areasToRemove, so
%   subdivisions are removed as well (e.g. 'LGd' also removes LGd-sh).

Coordinates(startsWith(Coordinates{:, 2}, areasToRemove), :) = [];
end
