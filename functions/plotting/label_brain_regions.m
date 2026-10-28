function label_brain_regions(region, y)
%LABEL_BRAIN_REGIONS Region name centred above each group of areas.
%   LABEL_BRAIN_REGIONS(region, y) writes, in the current axes at height y,
%   the name of each brain region over its areas (areas at x = 1, 3, 5, ...;
%   areas must be sorted by region).

st = brain_region_style();
regions = unique(region);
lengths = arrayfun(@(g) sum(region == g), regions)';
for i = 1:numel(regions)
    g = regions(i);
    text(sum(lengths(1:i)) * 2 - lengths(i), y, st.names(g), 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'Color', st.textColor(g, :));
end
end
