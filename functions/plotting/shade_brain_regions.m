function shade_brain_regions(region, yMax)
%SHADE_BRAIN_REGIONS Background band behind each area, coloured by brain region.
%   SHADE_BRAIN_REGIONS(region, yMax) draws, in the current axes, one band
%   from 0 to yMax per area (areas at x = 1, 3, 5, ...), coloured by its
%   brain region (parent group, see brain_region_style).

st = brain_region_style();
for k = 1:numel(region)
    jj = 2 * k - 1;
    g = region(k);
    p = patch([(jj-1) (jj-1) (jj+1) (jj+1)], [0 yMax yMax 0], st.shadeColor(g, :), ...
        'FaceAlpha', st.shadeAlpha(g), 'EdgeColor', 'none');
    hold on;
    uistack(p, 'bottom');
end
end
