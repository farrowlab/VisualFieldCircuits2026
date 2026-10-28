function E = experiment_list(celltype, expSet)
%EXPERIMENT_LIST Tracing experiments and injection site for each cell type.
%   E = EXPERIMENT_LIST(celltype, expSet) returns an N-by-2 cell array
%   {experiment name, 'Med' | 'Lat'}, where Med = medial sSC (upper visual
%   field) and Lat = lateral sSC (lower visual field) injection.
%
%   celltype: 'ntsr' or 'grp'
%   expSet:   'all'                     all brains
%             'noMidlineContamination'  (ntsr) brains without much labelling across the midline
%             'ipsiOnly'                (ntsr) brains whose injection-site labelling is
%                                       confined to the ipsilateral sSC
%             'retinal'                 (ntsr) brains used for the retinal contribution

switch celltype
    case 'ntsr'
        switch expSet
            case 'all'
                E = {'01660_ntsr_med',     'Med'
                     '01660_ntsr_lat',     'Lat'
                     '01662_ntsr_med',     'Med'
                     '01662_ntsr_lat',     'Lat'
                     '01830_ntsr_lat',     'Lat'
                     '01799_lp_ntsr_med',  'Med'
                     '01799_lp_ntsr_lat',  'Lat'
                     '01800_lp_ntsr_med',  'Med'
                     '01800_lp_ntsr_lat',  'Lat'
                     '01801_lp_ntsr_med',  'Med'
                     '01801_lp_ntsr_lat',  'Lat'
                     '01831_lp_ntsr_med',  'Med'
                     '01831_lp_ntsr_lat',  'Lat'
                     '01795_m02_ntsr_med', 'Med'};

            case 'noMidlineContamination'
                E = {'01660_ntsr_med',     'Med'
                     '01660_ntsr_lat',     'Lat'
                     '01662_ntsr_med',     'Med'
                     '01662_ntsr_lat',     'Lat'
                     '01830_ntsr_lat',     'Lat'
                     '01799_lp_ntsr_med',  'Med'
                     '01800_lp_ntsr_med',  'Med'
                     '01801_lp_ntsr_med',  'Med'
                     '01801_lp_ntsr_lat',  'Lat'
                     '01831_lp_ntsr_med',  'Med'
                     '01831_lp_ntsr_lat',  'Lat'
                     '01795_m02_ntsr_med', 'Med'};

            case 'ipsiOnly'
                E = {'01659_ntsr_med',     'Med'
                     '01660_ntsr_med',     'Med'
                     '01660_ntsr_lat',     'Lat'
                     '01662_ntsr_med',     'Med'
                     '01662_ntsr_lat',     'Lat'
                     '01799_lp_ntsr_lat',  'Lat'
                     '01800_lp_ntsr_med',  'Med'
                     '01800_lp_ntsr_lat',  'Lat'
                     '01801_lp_ntsr_lat',  'Lat'
                     '01831_lp_ntsr_lat',  'Lat'
                     '01795_m02_ntsr_med', 'Med'
                     '01830_ntsr_lat',     'Lat'};

            case 'retinal'
                E = {'01659_ntsr_med',     'Med'
                     '01659_ntsr_lat',     'Lat'
                     '01831_lp_ntsr_med',  'Med'
                     '01831_lp_ntsr_lat',  'Lat'};

            otherwise
                error('experiment_list:unknownSet', ...
                    'Unknown experiment set ''%s'' for %s.', expSet, celltype);
        end

    case 'grp'
        switch expSet
            case 'all'
                E = {'01976_grp_med',     'Med'
                     '01976_grp_lat',     'Lat'
                     '01977_grp_med',     'Med'
                     '01977_grp_lat',     'Lat'
                     '02017_grp_med',     'Med'
                     '02269_m02_grp_med', 'Med'
                     '02269_m02_grp_lat', 'Lat'
                     '02269_m03_grp_med', 'Med'
                     '02269_m03_grp_lat', 'Lat'
                     '02269_m04_grp_lat', 'Lat'
                     '02269_m05_grp_med', 'Med'};
            otherwise
                error('experiment_list:unknownSet', ...
                    'Unknown experiment set ''%s'' for %s.', expSet, celltype);
        end

    otherwise
        error('experiment_list:unknownCelltype', 'Unknown cell type ''%s''.', celltype);
end
end
