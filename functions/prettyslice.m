function prettyslice(slice,st,ax)

plotBrainOutlinesByAxis(slice,gca, [0.5 0.5 0.5],2.5); %add brain outline behind everything
%add brain area color
% brain_colors={[0.9568, 0.9294, 0.9019]};
brain_colors={[1, 0.9647 0.9333]};
sliceOutlineWithRegionVec(slice, 1, brain_colors, ax,[0.7 0.7 0.7], 0.5); hold on;

%add fiber tract color
fibers=uint16(st.sphinx_id(find(contains(st.structure_id_path, '1009'))));
fiber_colors=repmat({[0.75, 0.75, 0.75]},1,length(fibers));
sliceOutlineWithRegionVec(slice, fibers, fiber_colors, ax,'none', 0.5); hold on;

% %add supra-callosal cerebral white matter
% scwm=1190;
% scwm_colors=repmat({[0.75, 0.75, 0.75]},1,length(scwm));
% sliceOutlineWithRegionVec(slice, scwm, scwm_colors, ax,'none', 0.5); hold on;

%add medulla color
medulla=uint16(st.sphinx_id(find(contains(st.structure_id_path, '354'))));
medulla_colors=repmat({[0.75, 0.75, 0.75]},1,length(medulla));
sliceOutlineWithRegionVec(slice, medulla, medulla_colors, ax,'none', 0.5); hold on;

%add hippocampal formation color
hpf=455; %ventricle ids
hpf_colors=repmat({[0.75, 0.75, 0.75]},1,length(hpf));
sliceOutlineWithRegionVec(slice, hpf, hpf_colors, ax,'none', 0.5);

%add ventricle color
ventricle=1293:1304; %ventricle ids
ventricle_colors=repmat({[0.5, 0.5, 0.5]},1,length(ventricle));
sliceOutlineWithRegionVec(slice, ventricle, ventricle_colors, ax,'none', 0.5);

set(gca, 'YDir', 'reverse');
end