%% --- prepare data for heatmap of distribution overlap between med and lat NEW

clc;clear data;clear ttd;clear ttm; clear OD;
%load(fullfile(spath,'OD.mat')); load(fullfile(spath,'Coord.mat'));
celltypes = {'ntsr','grp'};
list_id = [];myexp_labels = {};nn = 0;n_Data=[];OD={}; Coord={};tt=0;
n_total={};
for c=1:length(celltypes)
    
    OD{c}=load(fullfile(spath,['OD_',celltypes{c},'.mat']));
    Coord{c}=load(fullfile(spath,['Coord_',celltypes{c},'.mat']));
    
    myType = flip(unique(OD{c}.OD_cell.type),2);
        
    for t = 1: length(myType)
        typeName = myType{t};   
        tableid = find(contains(OD{c}.OD_cell.type,typeName));
        list_id = [list_id tableid];
        tt=tt+1;
        data(tt).type = typeName;
        data(tt).table = load(fullfile(spath,[celltypes{c}, '_', typeName,'.mat']));
        thr = []; thr = data(tt).table.Threshold;
        n_Data(c,t)=size(data(tt).table.CoordMerge,2); 
        for j = 1:length(tableid)
            nn = nn+1;
            myexp_labels{nn,1} = [typeName,OD{c}.OD_cell.group{tableid(j)}];
            myexp_labels{nn,2} = OD{c}.OD_cell.exp{tableid(j)};
            n_total{c,t}(j)=OD{c}.OD_cell.nTotalminusSC{tableid(j)};
        end
    end
end
set(0,'defaultAxesFontSize',10);

n_Data=reshape(n_Data.', 1, []);
totalMed={}; totalLat={};
countLat={};densityLat={};percentLat={};threshlLat={};errPerLat={};errDenLat={};
countMed={};densityMed={};percentMed={};threshlMed={};errPerMed={};errDenMed={};
ll = length(BBB);
labels={};
% read Lat & Med data
for t = 1:4 %2 if 1 celltype, 4 if 2
    labels{t} = BBB(:,[1,3,6,8]);
    density = zeros(ll,n_Data(1,t));percent = zeros(ll,n_Data(1,t));
    dataCur = data(t).table.TT;data_tmp = [];
    threCur = data(t).table.Threshold;
    percent = zeros(ll,n_Data(1,t)); density = zeros(ll,n_Data(1,t));
    nrcount = zeros(ll,n_Data(1,t));
    for i = 1:ll
        data_tmp = [cell2mat(dataCur(i,7))];
        if ~isempty(data_tmp)
            nrcount(i,:) = data_tmp;
            percent(i,:) = [cell2mat(dataCur(i,9))];%percentage based on total inputs - SCs
            density(i,:) = [cell2mat(dataCur(i,15))];%density normalized to sum of all densitiess exclusing SCs
        end
    end

    % density_err = std(density,[],2);%standard deviation
    density_err = std(density,[],2)./sqrt(n_Data);%standard error of the mean (SEM)
    percent_err = std(percent,[],2)./sqrt(n_Data);%standard error of the mean (SEM)
    % percent_err = std(percent,[],2);%standard deviation

    % data
    if mod(t, 2) == 0       
        tt=t/2;
        countLat{tt} = nrcount;densityLat{tt} = density;
        percentLat{tt} = percent;threshlLat{tt} = threCur;
        errPerLat{tt} = percent_err;
        errDenLat{tt} = density_err;
    else
        tt=1+0.5*(t-1);
        countMed{tt} = nrcount;densityMed{tt} = density;
        percentMed{tt} = percent;threshlMed{tt} = threCur;
        errPerMed{tt} = percent_err;
        errDenMed{tt} = density_err;
    end

end

% remove SC
idSC = find(startsWith(labels{1,1}(:,1),'SCs'));%removes all SC layers from analysis
%remove MBsen,MBmot,MBsta
idMBdiv=find(ismember(labels{1,1}(:,1),["MBsen", "MBmot", "MBsta"]));
%filter out ctx
idnonCTX = find(~(cell2mat(labels{1,1}(:,3))==1));
%filter out IC subdivisions
idICdiv=find(ismember(labels{1,1}(:,1),["ICc", "ICd", "ICe", "ACAd", "VISam"]));


% --- Unified filtering: thresholds combined across celltypes, then applied to all ---
dmLat={}; dmMed= {}; pmLat={}; pmMed={}; 

% Thresholds Ntsr
thresT1 = 0;    % threshl condition
thresT2_ct1 = 1.6; % density condition
thresT3_ct1 = 1.11; % percent condition

% % Thresholds Grp
thresT1 = 0;    % threshl condition
thresT2_ct2 = 1.6; % density condition
thresT3_ct2 = 1.11; % percent condition

% general threshold
thresT1 = 0;    % threshl condition
thresT2 = 1.6; % density condition
thresT3 = 1.11; % percent condition


nCellTypes = length(celltypes);
nAreas0 = size(labels{1},1);

% Compute mean density and percent per area per celltype
for c = 1:nCellTypes
    dmLat{c} = mean(densityLat{c},2);
    dmMed{c} = mean(densityMed{c},2);
    pmLat{c} = mean(percentLat{c},2);
    pmMed{c} = mean(percentMed{c},2);
end

% --- Step 1: build global masks (union across celltypes) ---
mask_threshl = false(nAreas0,1);
mask_dm      = false(nAreas0,1);
mask_pm      = false(nAreas0,1);

for c = 1:nCellTypes
    mask_threshl = mask_threshl | (threshlLat{c} ~= thresT1);
    mask_threshl = mask_threshl | (threshlMed{c} ~= thresT1);
    mask_dm      = mask_dm      | (dmLat{c} >= thresT2);
    mask_dm      = mask_dm      | (dmMed{c} >= thresT2);
    mask_pm      = mask_pm      | (pmLat{c} >= thresT3);
    mask_pm      = mask_pm      | (pmMed{c} >= thresT3);
end

%separate for each celltypes
mask_threshl_ct1 = false(nAreas0,1);
mask_dm_ct1      = false(nAreas0,1);
mask_pm_ct1      = false(nAreas0,1);

mask_threshl_ct2 = false(nAreas0,1);
mask_dm_ct2      = false(nAreas0,1);
mask_pm_ct2      = false(nAreas0,1);

%ntsr
mask_threshl_ct1 = mask_threshl_ct1 | (threshlLat{1} ~= thresT1);
mask_threshl_ct1 = mask_threshl_ct1 | (threshlMed{1} ~= thresT1);
mask_dm_ct1      = mask_dm_ct1      | (dmLat{1} >= thresT2_ct1);
mask_dm_ct1      = mask_dm_ct1      | (dmMed{1} >= thresT2_ct1);
mask_pm_ct1      = mask_pm_ct1      | (pmLat{1} >= thresT3_ct1);
mask_pm_ct1      = mask_pm_ct1      | (pmMed{1} >= thresT3_ct1);

%grp
mask_threshl_ct2 = mask_threshl_ct2 | (threshlLat{2} ~= thresT1);
mask_threshl_ct2 = mask_threshl_ct2 | (threshlMed{2} ~= thresT1);
mask_dm_ct2      = mask_dm_ct2      | (dmLat{2} >= thresT2_ct2);
mask_dm_ct2      = mask_dm_ct2      | (dmMed{2} >= thresT2_ct2);
mask_pm_ct2      = mask_pm_ct2      | (pmLat{2} >= thresT3_ct2);
mask_pm_ct2      = mask_pm_ct2      | (pmMed{2} >= thresT3_ct2);


% --- Step 2: combine all three global conditions ---
mask_keep = mask_threshl & mask_dm & mask_pm;
mask_keep_ct1= mask_threshl_ct1 & mask_dm_ct1 & mask_pm_ct1;
mask_keep_ct2= mask_threshl_ct2 & mask_dm_ct2 & mask_pm_ct2;

% Remove global excludes
mask_keep(idSC)    = false;
mask_keep(idMBdiv) = false;
mask_keep(idICdiv) = false;

mask_keep_ct1(idSC)    = false;
mask_keep_ct1(idMBdiv) = false;
mask_keep_ct1(idICdiv) = false;

mask_keep_ct2(idSC)    = false;
mask_keep_ct2(idMBdiv) = false;
mask_keep_ct2(idICdiv) = false;


% Indices to keep
keepIdx = find(mask_keep);
keepIdx_ct1 = find(mask_keep_ct1);
keepIdx_ct2 = find(mask_keep_ct2);

fprintf('Final kept areas: %d of %d\n', numel(keepIdx), nAreas0);
fprintf('Final kept areas: %d of %d\n', numel(keepIdx_ct1), nAreas0);
fprintf('Final kept areas: %d of %d\n', numel(keepIdx_ct2), nAreas0);
% 
% % --- Step 3: apply SAME filter to all celltypes ---
% for c = 1:nCellTypes
%     labels{c}    = labels{c}(keepIdx,:);
%     countLat{c}    = countLat{c}(keepIdx,:);
%     countMed{c}    = countMed{c}(keepIdx,:);
%     percentLat{c}  = percentLat{c}(keepIdx,:);
%     percentMed{c}  = percentMed{c}(keepIdx,:);
%     pmLat{c}        = pmLat{c}(keepIdx);
%     pmMed{c}        = pmMed{c}(keepIdx);
%     densityLat{c} = densityLat{c}(keepIdx,:);
%     densityMed{c} = densityMed{c}(keepIdx,:);
%     dmLat{c}        = dmLat{c}(keepIdx);
%     dmMed{c}        = dmMed{c}(keepIdx);
%     errDenLat{c}    = errDenLat{c}(keepIdx);
%     errDenMed{c}    = errDenMed{c}(keepIdx);
%     errPerLat{c}    = errPerLat{c}(keepIdx);
%     errPerMed{c}    = errPerMed{c}(keepIdx);
% 
%     disp(['Total number of included ', celltypes{c}, ' medial input cells: ' ,num2str(sum(countMed{c},'all'))]);
%     disp(['Total number of included ', celltypes{c}, ' lateral input cells: ' ,num2str(sum(countLat{c},'all'))]);
%     disp(['Average cell counts per ', celltypes{c}, ' brain (medial): ', num2str(round(mean(sum(countMed{c},1))))]);
%     disp(['Average cell counts per ', celltypes{c}, ' brain (lateral): ', num2str(round(mean(sum(countLat{c},1))))]);
% end

%ntsr
labels{1}    = labels{1}(keepIdx_ct1,:);
countLat{1}    = countLat{1}(keepIdx_ct1,:);
countMed{1}    = countMed{1}(keepIdx_ct1,:);
percentLat{1}  = percentLat{1}(keepIdx_ct1,:);
percentMed{1}  = percentMed{1}(keepIdx_ct1,:);
pmLat{1}        = pmLat{1}(keepIdx_ct1);
pmMed{1}        = pmMed{1}(keepIdx_ct1);
densityLat{1} = densityLat{1}(keepIdx_ct1,:);
densityMed{1} = densityMed{1}(keepIdx_ct1,:);
dmLat{1}        = dmLat{1}(keepIdx_ct1);
dmMed{1}        = dmMed{1}(keepIdx_ct1);
errDenLat{1}    = errDenLat{1}(keepIdx_ct1);
errDenMed{1}    = errDenMed{1}(keepIdx_ct1);
errPerLat{1}    = errPerLat{1}(keepIdx_ct1);
errPerMed{1}    = errPerMed{1}(keepIdx_ct1);

disp(['Total number of included ', celltypes{1}, ' medial input cells: ' ,num2str(sum(countMed{1},'all'))]);
disp(['Total number of included ', celltypes{1}, ' lateral input cells: ' ,num2str(sum(countLat{1},'all'))]);
disp(['Average cell counts per ', celltypes{1}, ' brain (medial): ', num2str(round(mean(sum(countMed{1},1))))]);
disp(['Average cell counts per ', celltypes{1}, ' brain (lateral): ', num2str(round(mean(sum(countLat{1},1))))]);

%grp
labels{2}    = labels{2}(keepIdx_ct2,:);
countLat{2}    = countLat{2}(keepIdx_ct2,:);
countMed{2}    = countMed{2}(keepIdx_ct2,:);
percentLat{2}  = percentLat{2}(keepIdx_ct2,:);
percentMed{2}  = percentMed{2}(keepIdx_ct2,:);
pmLat{2}        = pmLat{2}(keepIdx_ct2);
pmMed{2}        = pmMed{2}(keepIdx_ct2);
densityLat{2} = densityLat{2}(keepIdx_ct2,:);
densityMed{2} = densityMed{2}(keepIdx_ct2,:);
dmLat{2}        = dmLat{2}(keepIdx_ct2);
dmMed{2}        = dmMed{2}(keepIdx_ct2);
errDenLat{2}    = errDenLat{2}(keepIdx_ct2);
errDenMed{2}    = errDenMed{2}(keepIdx_ct2);
errPerLat{2}    = errPerLat{2}(keepIdx_ct2);
errPerMed{2}    = errPerMed{2}(keepIdx_ct2);

disp(['Total number of included ', celltypes{2}, ' medial input cells: ' ,num2str(sum(countMed{2},'all'))]);
disp(['Total number of included ', celltypes{2}, ' lateral input cells: ' ,num2str(sum(countLat{2},'all'))]);
disp(['Average cell counts per ', celltypes{2}, ' brain (medial): ', num2str(round(mean(sum(countMed{2},1))))]);
disp(['Average cell counts per ', celltypes{2}, ' brain (lateral): ', num2str(round(mean(sum(countLat{2},1))))]);


% all_labels=labels{1};
all_labels=[];
for c=1:length(celltypes)
    all_labels=[all_labels;labels{c}];
end
[~, ia]=unique(all_labels(:,1), 'stable');
union_labels=all_labels(ia,:);
n_include = size(union_labels,1);
% disp(['total areas = ',mat2str(ll),' - ',mat2str(length(unique_nanIDX))]);
disp(['Include ' ,mat2str(n_include),' areas for all circuits']);


%% calculate distribution overlaps of areas along three main brain axes

%distrubition ML=3; AP=1; DV=2;


areas = string(union_labels(:,1));

% Store data for both celltypes
dt = {};
for c = 1:length(celltypes)
    for t = 1:length(myType)
        typeName = myType{t};
        dt{c,t} = load(fullfile(spath, [celltypes{c}, '_', typeName, '.mat']));
    end
end

labels = dt{1,1}.TT(:,1:3); % Assuming labels are the same across celltypes

% Match coordinates for each celltype separately
myCoord = {};
for c = 1:length(celltypes)
    min_exp_nr = min(length(dt{c,1}.expid), length(dt{c,2}.expid));
    
    if length(dt{c,1}.expid) == min_exp_nr
        for i = 1:min_exp_nr
            id = dt{c,1}.expid{i}(1:9);
            index_2 = find(contains(dt{c,2}.expid, id));
            if index_2
                myCoord{c,1}(:,i) = dt{c,1}.CoordMerge(:, i);
                myCoord{c,2}(:,i) = dt{c,2}.CoordMerge(:, index_2);
            end
        end
    else
        for i = 1:min_exp_nr
            id = dt{c,2}.expid{i}(1:9);
            index_1 = find(contains(dt{c,1}.expid, id));
            if index_1
                myCoord{c,1}(:,i) = dt{c,1}.CoordMerge(:, index_1);
                myCoord{c,2}(:,i) = dt{c,2}.CoordMerge(:, i);
            end
        end
    end
end

% Extract coordinates of individual areas for each celltype
coord_fin = {};
for c = 1:length(celltypes)
    for a = 1:length(areas)
        coord_tot = {};
        for t = 1:length(myType)
            coll = col(t+2,:);
            idSC = find(ismember(labels(:,1), areas(a)));
            myCoord_temp = myCoord{c,t}(idSC,:);
            
            [sa, sl] = size(myCoord_temp);
            Coord_type = {};
            
            for i = 1:sl
                coordXYZ = cell2mat(myCoord_temp(:,i));
                if ~isempty(coordXYZ)
                    coordXYZ = coordXYZ(coordXYZ(:,3) >= 0, :); % keep only ipsilateral
                    coordXYZ = coordXYZ + randn(size(coordXYZ))/300; % add noise
                    coordXYZ(:,1) = coordXYZ(:,1) + (-0.05 + (0.05+0.05).*rand(size(coordXYZ,1),1)); % jitter AP
                    Coord_type{i} = coordXYZ;
                else
                    Coord_type{i} = [];
                end
            end
            coord_tot{1,t} = Coord_type;
        end
        coord_fin{c,a} = coord_tot; % Now indexed by celltype and area
    end
end

%% stability checks and significance
% Modified section of your code with all three axes
ML_ksd = []; ML_pval_perm = []; ML_stability = [];
AP_ksd = []; AP_pval_perm = []; AP_stability = [];
DV_ksd = []; DV_pval_perm = []; DV_stability = [];

min_cells_strict = 30;  % More conservative threshold
min_cells_permutation = 15;  % Minimum for permutation test

for a = 1:length(coord_fin)   
    for i = 1:length(coord_fin{1,1}{1,1})
        n_cells_type1 = size(coord_fin{1,a}{1,1}{i}, 1);
        n_cells_type2 = size(coord_fin{1,a}{1,2}{i}, 1);
        
        if n_cells_type1 >= min_cells_strict && n_cells_type2 >= min_cells_strict
            % Use KDE with stability check for all three axes
            
            %---ML axis (column 3)---
            [overlap_ML, p_ML, ~] = permutation_overlap_test(...
                coord_fin{1,a}{1,1}{i}, coord_fin{1,a}{1,2}{i}, 3, 1000);
            ML_ksd(a,i) = overlap_ML;
            ML_pval_perm(a,i) = 1-p_ML;
            
            % Assess ML stability
            stab_ML = assess_kde_stability(...
                coord_fin{1,a}{1,1}{i}(:,3), coord_fin{1,a}{1,2}{i}(:,3), 100);
            ML_stability(a,i) = stab_ML.cv;
            
            %---AP axis (column 1)---
            [overlap_AP, p_AP, ~] = permutation_overlap_test(...
                coord_fin{1,a}{1,1}{i}, coord_fin{1,a}{1,2}{i}, 1, 1000);
            AP_ksd(a,i) = overlap_AP;
            AP_pval_perm(a,i) = 1-p_AP;
            
            % Assess AP stability
            stab_AP = assess_kde_stability(...
                coord_fin{1,a}{1,1}{i}(:,1), coord_fin{1,a}{1,2}{i}(:,1), 100);
            AP_stability(a,i) = stab_AP.cv;
            
            %---DV axis (column 2)---
            [overlap_DV, p_DV, ~] = permutation_overlap_test(...
                coord_fin{1,a}{1,1}{i}, coord_fin{1,a}{1,2}{i}, 2, 1000);
            DV_ksd(a,i) = overlap_DV;
            DV_pval_perm(a,i) = 1-p_DV;
            
            % Assess DV stability
            stab_DV = assess_kde_stability(...
                coord_fin{1,a}{1,1}{i}(:,2), coord_fin{1,a}{1,2}{i}(:,2), 100);
            DV_stability(a,i) = stab_DV.cv;
            
        elseif n_cells_type1 >= min_cells_permutation && n_cells_type2 >= min_cells_permutation
            % Use permutation test only (no KDE overlap) for smaller samples
            
            % ML axis
            [~, p_ML] = kstest2(coord_fin{1,a}{1,1}{i}(:,3), ...
                               coord_fin{1,a}{1,2}{i}(:,3));
            ML_ksd(a,i) = NaN;  % Don't trust KDE overlap
            ML_pval_perm(a,i) = 1-p_ML;
            ML_stability(a,i) = NaN;
            
            % AP axis
            [~, p_AP] = kstest2(coord_fin{1,a}{1,1}{i}(:,1), ...
                               coord_fin{1,a}{1,2}{i}(:,1));
            AP_ksd(a,i) = NaN;
            AP_pval_perm(a,i) = 1-p_AP;
            AP_stability(a,i) = NaN;
            
            % DV axis
            [~, p_DV] = kstest2(coord_fin{1,a}{1,1}{i}(:,2), ...
                               coord_fin{1,a}{1,2}{i}(:,2));
            DV_ksd(a,i) = NaN;
            DV_pval_perm(a,i) = 1-p_DV;
            DV_stability(a,i) = NaN;
            
        else
            % Too few cells - skip analysis for all axes
            ML_ksd(a,i) = NaN;
            ML_pval_perm(a,i) = NaN;
            ML_stability(a,i) = NaN;
            
            AP_ksd(a,i) = NaN;
            AP_pval_perm(a,i) = NaN;
            AP_stability(a,i) = NaN;
            
            DV_ksd(a,i) = NaN;
            DV_pval_perm(a,i) = NaN;
            DV_stability(a,i) = NaN;
        end
    end
end

% Complete stability check for all three axes
fprintf('=== STABILITY ASSESSMENT ===\n\n');

% ML Stability
stability_ML = ML_stability(~isnan(ML_stability));
fprintf('ML Stability Assessment:\n');
fprintf('Mean CV: %.3f\n', mean(stability_ML));
fprintf('Areas with CV > 0.2 (unstable): %d/%d\n', ...
    sum(stability_ML > 0.2), length(stability_ML));
fprintf('Areas with CV > 0.3 (very unstable): %d/%d\n\n', ...
    sum(stability_ML > 0.3), length(stability_ML));

% AP Stability
stability_AP = AP_stability(~isnan(AP_stability));
fprintf('AP Stability Assessment:\n');
fprintf('Mean CV: %.3f\n', mean(stability_AP));
fprintf('Areas with CV > 0.2 (unstable): %d/%d\n', ...
    sum(stability_AP > 0.2), length(stability_AP));
fprintf('Areas with CV > 0.3 (very unstable): %d/%d\n\n', ...
    sum(stability_AP > 0.3), length(stability_AP));

% DV Stability
stability_DV = DV_stability(~isnan(DV_stability));
fprintf('DV Stability Assessment:\n');
fprintf('Mean CV: %.3f\n', mean(stability_DV));
fprintf('Areas with CV > 0.2 (unstable): %d/%d\n', ...
    sum(stability_DV > 0.2), length(stability_DV));
fprintf('Areas with CV > 0.3 (very unstable): %d/%d\n\n', ...
    sum(stability_DV > 0.3), length(stability_DV));

% Identify problematic areas across all axes
problematic_ML = areas(any(ML_stability > 0.2, 2));
problematic_AP = areas(any(AP_stability > 0.2, 2));
problematic_DV = areas(any(DV_stability > 0.2, 2));

if ~isempty(problematic_ML)
    fprintf('Areas with unstable ML KDE estimates:\n');
    disp(problematic_ML);
end
if ~isempty(problematic_AP)
    fprintf('Areas with unstable AP KDE estimates:\n');
    disp(problematic_AP);
end
if ~isempty(problematic_DV)
    fprintf('Areas with unstable DV KDE estimates:\n');
    disp(problematic_DV);
end

% Comprehensive decision matrix for all axes
fprintf('\n=== DECISION MATRIX ===\n');
decision_matrix = table();

for a = 1:length(areas)
    % Get average stability across animals for each axis
    ml_stab = mean(ML_stability(a,:), 'omitnan');
    ap_stab = mean(AP_stability(a,:), 'omitnan');
    dv_stab = mean(DV_stability(a,:), 'omitnan');
    
    % Get average sample sizes for both cell types
    avg_n1 = 0; avg_n2 = 0;
    count = 0;
    for i = 1:size(coord_fin{1,1}{1,1}, 2)
        if ~isempty(coord_fin{1,a}{1,1}{i}) && ~isempty(coord_fin{1,a}{1,2}{i})
            avg_n1 = avg_n1 + size(coord_fin{1,a}{1,1}{i}, 1);
            avg_n2 = avg_n2 + size(coord_fin{1,a}{1,2}{i}, 1);
            count = count + 1;
        end
    end
    if count > 0
        avg_n1 = avg_n1 / count;
        avg_n2 = avg_n2 / count;
        min_avg_n = min(avg_n1, avg_n2);
    else
        min_avg_n = 0;
    end
    
    % Make recommendation based on worst-case stability
    max_cv = max([ml_stab, ap_stab, dv_stab]);
    if min_avg_n >= 30 && max_cv < 0.2
        recommendation = "Good - Use KDE";
    elseif min_avg_n >= 20 && max_cv < 0.3
        recommendation = "Marginal - Consider alternatives";
    elseif min_avg_n >= 15
        recommendation = "Poor - Use permutation test only";
    else
        recommendation = "Insufficient data";
    end
    
    decision_matrix = [decision_matrix; 
        table(areas(a), round(min_avg_n), ...
              round(ml_stab,3), round(ap_stab,3), round(dv_stab,3), ...
              round(max_cv,3), recommendation)];
end

decision_matrix.Properties.VariableNames = ...
    {'Area', 'Min_Avg_N', 'ML_CV', 'AP_CV', 'DV_CV', 'Max_CV', 'Recommendation'};
disp(decision_matrix);

% Summary statistics table
fprintf('\n=== SUMMARY STATISTICS ===\n');

% Clean the data as in your original code
mask = isnan(ML_ksd) | (ML_ksd==0);
ML_ksd_clean = ML_ksd;
AP_ksd_clean = AP_ksd;
DV_ksd_clean = DV_ksd;
ML_ksd_clean(mask) = nan;
AP_ksd_clean(mask) = nan;
DV_ksd_clean(mask) = nan;

% Calculate means
ML_mean = mean(ML_ksd_clean, 2, 'omitnan');
AP_mean = mean(AP_ksd_clean, 2, 'omitnan');
DV_mean = mean(DV_ksd_clean, 2, 'omitnan');

% Calculate mean p-values
ML_pval_mean = mean(ML_pval_perm, 2, 'omitnan');
AP_pval_mean = mean(AP_pval_perm, 2, 'omitnan');
DV_pval_mean = mean(DV_pval_perm, 2, 'omitnan');

% Create final results table
final_results = table(areas, ...
    round(ML_mean, 3), round(AP_mean, 3), round(DV_mean, 3), ...
    round(ML_pval_mean, 4), round(AP_pval_mean, 4), round(DV_pval_mean, 4));
final_results.Properties.VariableNames = ...
    {'Area', 'ML_Overlap', 'AP_Overlap', 'DV_Overlap', ...
     'ML_pval', 'AP_pval', 'DV_pval'};

fprintf('\nFinal overlap results with p-values:\n');
disp(final_results);

% Count significant results after multiple comparison correction
alpha = 0.05;
fprintf('\n=== PERMUTATION TEST RESULTS ===\n');
fprintf('Significant differences (p < %.3f):\n', alpha);
fprintf('ML: %d/%d areas (%.1f%%)\n', ...
    sum(ML_pval_mean < alpha), sum(~isnan(ML_pval_mean)), ...
    100*sum(ML_pval_mean < alpha)/sum(~isnan(ML_pval_mean)));
fprintf('AP: %d/%d areas (%.1f%%)\n', ...
    sum(AP_pval_mean < alpha), sum(~isnan(AP_pval_mean)), ...
    100*sum(AP_pval_mean < alpha)/sum(~isnan(AP_pval_mean)));
fprintf('DV: %d/%d areas (%.1f%%)\n', ...
    sum(DV_pval_mean < alpha), sum(~isnan(DV_pval_mean)), ...
    100*sum(DV_pval_mean < alpha)/sum(~isnan(DV_pval_mean)));


% Bootstrap stability test for your KDE overlap measures
function stability = assess_kde_stability(coord1, coord2, n_bootstrap)
    if nargin < 3
        n_bootstrap = 100;
    end
    
    overlaps = zeros(n_bootstrap, 1);
    n1 = size(coord1, 1);
    n2 = size(coord2, 1);
    
    for b = 1:n_bootstrap
        % Resample with replacement
        idx1 = randsample(n1, n1, true);
        idx2 = randsample(n2, n2, true);
        
        % Calculate overlap with bootstrap samples
        [y1, ~] = ksdensity(coord1(idx1));
        x1 = linspace(min(coord1(idx1)), max(coord1(idx1)), 100);
        [y2, ~] = ksdensity(coord2(idx2));
        x2 = linspace(min(coord2(idx2)), max(coord2(idx2)), 100);
        
        % Common x-axis
        newx = linspace(min([min(coord1(idx1)), min(coord2(idx2))]), ...
                       max([max(coord1(idx1)), max(coord2(idx2))]), 100);
        new_y1 = interp1(x1, y1, newx);
        new_y2 = interp1(x2, y2, newx);
        
        pmin = min(new_y1, new_y2);
        pmax = max(new_y1, new_y2);
        overlaps(b) = sum(pmin) ./ sum(pmax);
    end
    
    stability.mean_overlap = mean(overlaps);
    stability.std_overlap = std(overlaps);
    stability.cv = std(overlaps) / mean(overlaps); % Coefficient of variation
    stability.ci = prctile(overlaps, [2.5, 97.5]); % 95% CI
end

function [overlap_obs, p_value, null_dist] = permutation_overlap_test(coord1, coord2, axis_idx, n_perm)
    if nargin < 4
        n_perm = 1000;
    end
    
    % Extract the axis of interest (1=AP, 2=DV, 3=ML)
    data1 = coord1(:, axis_idx);
    data2 = coord2(:, axis_idx);
    
    % Calculate observed overlap
    overlap_obs = calculate_overlap(data1, data2);
    
    % Combine all data
    all_data = [data1; data2];
    n1 = length(data1);
    n2 = length(data2);
    n_total = n1 + n2;
    
    % Permutation test
    null_dist = zeros(n_perm, 1);
    for p = 1:n_perm
        % Randomly permute combined data
        perm_idx = randperm(n_total);
        perm_data1 = all_data(perm_idx(1:n1));
        perm_data2 = all_data(perm_idx(n1+1:end));
        
        % Calculate overlap for permuted data
        null_dist(p) = calculate_overlap(perm_data1, perm_data2);
    end
    
    % Calculate p-value (two-tailed)
    % Tests if observed overlap is significantly different from random
    p_value = sum(abs(null_dist - 0.5) >= abs(overlap_obs - 0.5)) / n_perm;
end

function overlap = calculate_overlap(data1, data2)
    % Your existing overlap calculation
    [y1, ~] = ksdensity(data1);
    x1 = linspace(min(data1), max(data1), 100);
    [y2, ~] = ksdensity(data2);
    x2 = linspace(min(data2), max(data2), 100);
    
    newx = linspace(min([min(data1), min(data2)]), ...
                   max([max(data1), max(data2)]), 100);
    new_y1 = interp1(x1, y1, newx);
    new_y2 = interp1(x2, y2, newx);
    
    pmin = min(new_y1, new_y2);
    pmax = max(new_y1, new_y2);
    overlap = sum(pmin) ./ sum(pmax);
end



% Create a diagnostic plot
figure('Position', [100, 100, 1200, 400]);

% For each axis (ML, AP, DV)
axes_names = {'ML', 'AP', 'DV'};
stability_data = {ML_stability, AP_stability, DV_stability};

for ax = 1:3
    subplot(1, 3, ax);
    
    % Extract sample sizes and stability for plotting
    sample_sizes = [];
    stability_values = [];
    
    for a = 1:length(coord_fin)
        for i = 1:size(coord_fin{1,1}{1,1}, 2)
            if ~isnan(stability_data{ax}(a,i))
                n1 = size(coord_fin{1,a}{1,1}{i}, 1);
                n2 = size(coord_fin{1,a}{1,2}{i}, 1);
                min_n = min(n1, n2);
                sample_sizes = [sample_sizes, min_n];
                stability_values = [stability_values, stability_data{ax}(a,i)];
            end
        end
    end
    
    % Plot
    scatter(sample_sizes, stability_values, 50, 'filled');
    hold on;
    
    % Add reference lines
    yline(0.1, 'g--', 'CV=0.1 (Good)', 'LineWidth', 1.5);
    yline(0.2, 'y--', 'CV=0.2 (Acceptable)', 'LineWidth', 1.5);
    yline(0.3, 'r--', 'CV=0.3 (Poor)', 'LineWidth', 1.5);
    
    xlabel('Minimum Sample Size');
    ylabel('Coefficient of Variation');
    title(sprintf('%s Axis Stability', axes_names{ax}));
    grid on;
    
    % Add vertical line at n=20 and n=30
    xline(20, 'k:', 'n=20');
    xline(30, 'k:', 'n=30');
end
subtitle('KDE Stability Assessment');

% Examine significant differences after permutation test
alpha = 0.05;

% Count significant results
sig_ML = sum(ML_pval_perm < alpha, 'all', 'omitnan');
sig_AP = sum(AP_pval_perm < alpha, 'all', 'omitnan');
sig_DV = sum(DV_pval_perm < alpha, 'all', 'omitnan');

total_ML = sum(~isnan(ML_pval_perm), 'all');
total_AP = sum(~isnan(AP_pval_perm), 'all');
total_DV = sum(~isnan(DV_pval_perm), 'all');

fprintf('\nPermutation Test Results (p < 0.05):\n');
fprintf('ML: %d/%d significant (%.1f%%)\n', sig_ML, total_ML, 100*sig_ML/total_ML);
fprintf('AP: %d/%d significant (%.1f%%)\n', sig_AP, total_AP, 100*sig_AP/total_AP);
fprintf('DV: %d/%d significant (%.1f%%)\n', sig_DV, total_DV, 100*sig_DV/total_DV);

% Plot overlap values with confidence based on stability
figure('Position', [100, 100, 800, 600]);

% Combine all axes data
all_overlaps = [ML_ksd(:); AP_ksd(:); DV_ksd(:)];
all_stability = [ML_stability(:); AP_stability(:); DV_stability(:)];

% Remove NaNs
valid_idx = ~isnan(all_overlaps) & ~isnan(all_stability);
all_overlaps = all_overlaps(valid_idx);
all_stability = all_stability(valid_idx);

% Color code by stability
scatter(all_overlaps, 1:length(all_overlaps), 50, all_stability, 'filled');
colormap(flipud(hot));
c = colorbar;
c.Label.String = 'CV (lower is better)';
xlabel('Overlap Value');
ylabel('Measurement Index');
title('Distribution Overlaps Colored by Stability');

% Add vertical line at 0.5 (random expectation)
xline(0.5, 'k--', 'Random Overlap');


% Create a decision matrix for each area
decision_matrix = table();
for a = 1:length(areas)
    % Get average stability across animals
    ml_stab = mean(ML_stability(a,:), 'omitnan');
    ap_stab = mean(AP_stability(a,:), 'omitnan');
    dv_stab = mean(DV_stability(a,:), 'omitnan');
    
    % Get average sample size
    avg_n = 0;
    count = 0;
    for i = 1:size(coord_fin{1,1}{1,1}, 2)
        if ~isempty(coord_fin{1,a}{1,1}{i})
            avg_n = avg_n + size(coord_fin{1,a}{1,1}{i}, 1);
            count = count + 1;
        end
    end
    avg_n = avg_n / count;
    
    % Make recommendation
    if avg_n >= 30 && max([ml_stab, ap_stab, dv_stab]) < 0.2
        recommendation = "Good - Use KDE";
    elseif avg_n >= 20 && max([ml_stab, ap_stab, dv_stab]) < 0.3
        recommendation = "Marginal - Consider alternatives";
    else
        recommendation = "Poor - Need more data";
    end
    
    decision_matrix = [decision_matrix; 
        table(areas(a), avg_n, ml_stab, ap_stab, dv_stab, recommendation)];
end

decision_matrix.Properties.VariableNames = {'Area', 'Avg_N', 'ML_CV', 'AP_CV', 'DV_CV', 'Recommendation'};
disp(decision_matrix);

%% permutation analysis and significance
% Complete Permutation Analysis for All Areas and Directions
% This replaces the section after coord_fin is created

% Parameters
n_perm = 500;  % Number of permutations (increase to 1000 for publication)
min_cells = 20;  % Minimum cells required

% Initialize storage matrices
n_areas = length(areas);
max_animals = size(coord_fin{1}{1}, 2);

% Initialize results matrices
ML_overlap_all = NaN(n_areas, max_animals);
AP_overlap_all = NaN(n_areas, max_animals);
DV_overlap_all = NaN(n_areas, max_animals);

ML_pval_all = NaN(n_areas, max_animals);
AP_pval_all = NaN(n_areas, max_animals);
DV_pval_all = NaN(n_areas, max_animals);

% Run permutation tests for all areas and animals
fprintf('Running permutation tests for all areas and directions...\n');
fprintf('Progress: ');

for a = 1:n_areas
    fprintf('%d/%d...', a, n_areas);
    
    for i = 1:size(coord_fin{a}{1}, 2)
        % Check if data exists for both cell types
        if isempty(coord_fin{a}{1}{i}) || isempty(coord_fin{a}{2}{i})
            continue;
        end
        
        % Check sample sizes
        n1 = size(coord_fin{a}{1}{i}, 1);
        n2 = size(coord_fin{a}{2}{i}, 1);
        
        if n1 >= min_cells && n2 >= min_cells
            % ML axis (column 3)
            [ML_overlap_all(a,i), ML_pval_all(a,i)] = ...
                quick_permutation_test(coord_fin{a}{1}{i}(:,3), ...
                                     coord_fin{a}{2}{i}(:,3), n_perm);
            
            % AP axis (column 1)
            [AP_overlap_all(a,i), AP_pval_all(a,i)] = ...
                quick_permutation_test(coord_fin{a}{1}{i}(:,1), ...
                                     coord_fin{a}{2}{i}(:,1), n_perm);
            
            % DV axis (column 2)
            [DV_overlap_all(a,i), DV_pval_all(a,i)] = ...
                quick_permutation_test(coord_fin{a}{1}{i}(:,2), ...
                                     coord_fin{a}{2}{i}(:,2), n_perm);
        end
    end
end
fprintf(' Done!\n');

% Calculate mean values across animals for each area
ML_overlap_mean = mean(ML_overlap_all, 2, 'omitnan');
AP_overlap_mean = mean(AP_overlap_all, 2, 'omitnan');
DV_overlap_mean = mean(DV_overlap_all, 2, 'omitnan');

ML_pval_mean = mean(ML_pval_all, 2, 'omitnan');
AP_pval_mean = mean(AP_pval_all, 2, 'omitnan');
DV_pval_mean = mean(DV_pval_all, 2, 'omitnan');

% Count number of animals per area
n_animals_per_area = sum(~isnan(ML_overlap_all), 2);

% Create comprehensive visualization
create_all_plots(areas, myType, ML_overlap_mean, AP_overlap_mean, DV_overlap_mean, ...
                 ML_pval_mean, AP_pval_mean, DV_pval_mean, n_animals_per_area, ...
                 coord_fin, ML_overlap_all, AP_overlap_all, DV_overlap_all);

% Create results table
results_table = table(areas, n_animals_per_area, ...
    round(ML_overlap_mean, 3), round(AP_overlap_mean, 3), round(DV_overlap_mean, 3), ...
    round(ML_pval_mean, 4), round(AP_pval_mean, 4), round(DV_pval_mean, 4));

results_table.Properties.VariableNames = {'Area', 'N_Animals', ...
    'ML_Overlap', 'AP_Overlap', 'DV_Overlap', ...
    'ML_pval', 'AP_pval', 'DV_pval'};

% Display results
fprintf('\n=== PERMUTATION TEST RESULTS ===\n');
disp(results_table);

% Summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
alpha = 0.05;
fprintf('Significant differences (p < %.3f):\n', alpha);
fprintf('ML: %d/%d areas (%.1f%%)\n', ...
    sum(ML_pval_mean < alpha & ~isnan(ML_pval_mean)), sum(~isnan(ML_pval_mean)), ...
    100*sum(ML_pval_mean < alpha & ~isnan(ML_pval_mean))/sum(~isnan(ML_pval_mean)));
fprintf('AP: %d/%d areas (%.1f%%)\n', ...
    sum(AP_pval_mean < alpha & ~isnan(AP_pval_mean)), sum(~isnan(AP_pval_mean)), ...
    100*sum(AP_pval_mean < alpha & ~isnan(AP_pval_mean))/sum(~isnan(AP_pval_mean)));
fprintf('DV: %d/%d areas (%.1f%%)\n', ...
    sum(DV_pval_mean < alpha & ~isnan(DV_pval_mean)), sum(~isnan(DV_pval_mean)), ...
    100*sum(DV_pval_mean < alpha & ~isnan(DV_pval_mean))/sum(~isnan(DV_pval_mean)));

% Save results to file (optional)
save('permutation_test_results.mat', 'results_table', ...
     'ML_overlap_all', 'AP_overlap_all', 'DV_overlap_all', ...
     'ML_pval_all', 'AP_pval_all', 'DV_pval_all');

% FUNCTION DEFINITIONS

function [overlap_obs, p_value] = quick_permutation_test(data1, data2, n_perm)
    % Calculate observed overlap
    overlap_obs = calculate_kde_overlap(data1, data2);
    
    % Combine data for permutation
    all_data = [data1(:); data2(:)];
    n1 = length(data1);
    n_total = length(all_data);
    
    % Permutation test
    null_dist = zeros(n_perm, 1);
    for p = 1:n_perm
        % Shuffle and split
        perm_idx = randperm(n_total);
        perm_data1 = all_data(perm_idx(1:n1));
        perm_data2 = all_data(perm_idx(n1+1:end));
        
        % Calculate permuted overlap
        null_dist(p) = calculate_kde_overlap(perm_data1, perm_data2);
    end
    
    % Two-tailed p-value
    p_value = 1-(sum(abs(null_dist - 0.5) >= abs(overlap_obs - 0.5)) / n_perm);
end

function overlap = calculate_kde_overlap(data1, data2)
    try
        % Get densities
        [f1, x1] = ksdensity(data1);
        [f2, x2] = ksdensity(data2);
        
        % Common x-axis
        x_common = linspace(min([min(x1), min(x2)]), ...
                          max([max(x1), max(x2)]), 100);
        
        % Interpolate to common axis
        f1_interp = interp1(x1, f1, x_common, 'linear', 0);
        f2_interp = interp1(x2, f2, x_common, 'linear', 0);
        
        % Calculate overlap
        overlap = sum(min(f1_interp, f2_interp)) / sum(max(f1_interp, f2_interp));
    catch
        overlap = NaN;
    end
end

function create_all_plots(areas, myType, ML_overlap, AP_overlap, DV_overlap, ...
                         ML_pval, AP_pval, DV_pval, n_animals, coord_fin, ...
                         ML_all, AP_all, DV_all)
    
    % Figure 1: Main results
    figure('Position', [50, 50, 1500, 900], 'Name', 'Permutation Test Results - All Areas');
    
    n_areas = length(areas);
    x = 1:n_areas;
    
    % Remove NaN areas for plotting
    valid_idx = ~isnan(ML_overlap);
    x_valid = x(valid_idx);
    
    % Subplot 1: Overlap values for all axes
    subplot(2,3,1:3);
    hold on;
    
    % Plot overlaps
    h1 = plot(x_valid, ML_overlap(valid_idx), 'o-', 'LineWidth', 2, ...
              'MarkerSize', 8, 'Color', [0.2 0.4 0.8], 'DisplayName', 'ML');
    h2 = plot(x_valid, AP_overlap(valid_idx), 's-', 'LineWidth', 2, ...
              'MarkerSize', 8, 'Color', [0.8 0.3 0.3], 'DisplayName', 'AP');
    h3 = plot(x_valid, DV_overlap(valid_idx), '^-', 'LineWidth', 2, ...
              'MarkerSize', 8, 'Color', [0.3 0.7 0.3], 'DisplayName', 'DV');
    
    % Mark significant differences
    sig_ml = find(ML_pval < 0.05 & valid_idx);
    sig_ap = find(AP_pval < 0.05 & valid_idx);
    sig_dv = find(DV_pval < 0.05 & valid_idx);
    
    if ~isempty(sig_ml)
        plot(sig_ml, ML_overlap(sig_ml), 'ko', 'MarkerSize', 12, 'LineWidth', 2);
    end
    if ~isempty(sig_ap)
        plot(sig_ap, AP_overlap(sig_ap), 'ks', 'MarkerSize', 12, 'LineWidth', 2);
    end
    if ~isempty(sig_dv)
        plot(sig_dv, DV_overlap(sig_dv), 'k^', 'MarkerSize', 12, 'LineWidth', 2);
    end
    
    % Reference line
    yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    
    % Formatting
    xlabel('Brain Area', 'FontSize', 12);
    ylabel('Distribution Overlap', 'FontSize', 12);
    title(sprintf('Distribution Overlaps: %s vs %s (Black outline = p < 0.05)', ...
          myType{1}, myType{2}), 'FontSize', 14);
    xticks(x_valid);
    xticklabels(areas(valid_idx));
    xtickangle(45);
    legend([h1, h2, h3], 'Location', 'best', 'FontSize', 10);
    grid on;
    ylim([0 1]);
    xlim([0.5, n_areas+0.5]);
    
    % Subplot 2: P-values
    subplot(2,3,4);
    hold on;
    bar_data = [ML_pval(valid_idx), AP_pval(valid_idx), DV_pval(valid_idx)];
    b = bar(bar_data);
    b(1).FaceColor = [0.2 0.4 0.8];
    b(2).FaceColor = [0.8 0.3 0.3];
    b(3).FaceColor = [0.3 0.7 0.3];
    
    yline(0.05, 'r--', 'LineWidth', 2);
    xlabel('Brain Area', 'FontSize', 10);
    ylabel('p-value', 'FontSize', 10);
    title('Permutation Test P-values', 'FontSize', 12);
    xticks(1:sum(valid_idx));
    xticklabels(areas(valid_idx));
    xtickangle(45);
    legend({'ML', 'AP', 'DV', 'p=0.05'}, 'Location', 'best', 'FontSize', 9);
    grid on;
    
    % Subplot 3: Effect sizes
    subplot(2,3,5);
    effect_ml = abs(ML_overlap(valid_idx) - 0.5);
    effect_ap = abs(AP_overlap(valid_idx) - 0.5);
    effect_dv = abs(DV_overlap(valid_idx) - 0.5);
    
    bar_data = [effect_ml, effect_ap, effect_dv];
    b = bar(bar_data);
    b(1).FaceColor = [0.2 0.4 0.8];
    b(2).FaceColor = [0.8 0.3 0.3];
    b(3).FaceColor = [0.3 0.7 0.3];
    
    xlabel('Brain Area', 'FontSize', 10);
    ylabel('|Overlap - 0.5|', 'FontSize', 10);
    title('Effect Size', 'FontSize', 12);
    xticks(1:sum(valid_idx));
    xticklabels(areas(valid_idx));
    xtickangle(45);
    legend({'ML', 'AP', 'DV'}, 'Location', 'best', 'FontSize', 9);
    grid on;
    
    % Subplot 4: Sample sizes
    subplot(2,3,6);
    bar(n_animals(valid_idx), 'FaceColor', [0.5 0.5 0.5]);
    xlabel('Brain Area', 'FontSize', 10);
    ylabel('N Animals', 'FontSize', 10);
    title('Sample Sizes', 'FontSize', 12);
    xticks(1:sum(valid_idx));
    xticklabels(areas(valid_idx));
    xtickangle(45);
    grid on;
    
    % Figure 2: Heatmap view
    figure('Position', [100, 100, 1200, 800], 'Name', 'Overlap Heatmaps');
    
    % Create heatmaps for each axis
    subplot(2,2,1);
    imagesc(ML_all);
    colorbar;
    colormap(jet);
    caxis([0 1]);
    xlabel('Animal', 'FontSize', 10);
    ylabel('Brain Area', 'FontSize', 10);
    title('ML Axis Overlaps', 'FontSize', 12);
    yticks(1:n_areas);
    yticklabels(areas);
    
    subplot(2,2,2);
    imagesc(AP_all);
    colorbar;
    colormap(jet);
    caxis([0 1]);
    xlabel('Animal', 'FontSize', 10);
    ylabel('Brain Area', 'FontSize', 10);
    title('AP Axis Overlaps', 'FontSize', 12);
    yticks(1:n_areas);
    yticklabels(areas);
    
    subplot(2,2,3);
    imagesc(DV_all);
    colorbar;
    colormap(jet);
    caxis([0 1]);
    xlabel('Animal', 'FontSize', 10);
    ylabel('Brain Area', 'FontSize', 10);
    title('DV Axis Overlaps', 'FontSize', 12);
    yticks(1:n_areas);
    yticklabels(areas);
    
    % Summary in subplot 4
    subplot(2,2,4);
    axis off;
    
    % Calculate summary statistics
    alpha = 0.05;
    n_sig_ml = sum(ML_pval < alpha & ~isnan(ML_pval));
    n_sig_ap = sum(AP_pval < alpha & ~isnan(AP_pval));
    n_sig_dv = sum(DV_pval < alpha & ~isnan(DV_pval));
    n_total = sum(~isnan(ML_pval));
    
    % Find strongest effects
    [max_ml, idx_ml] = max(effect_ml);
    [max_ap, idx_ap] = max(effect_ap);
    [max_dv, idx_dv] = max(effect_dv);
    
    valid_areas = areas(valid_idx);
    
    summary_text = sprintf([...
        '\\fontsize{11}Summary Statistics\n\n' ...
        '\\fontsize{10}Significant Areas (p<0.05):\n' ...
        'ML: %d/%d (%.1f%%)\n' ...
        'AP: %d/%d (%.1f%%)\n' ...
        'DV: %d/%d (%.1f%%)\n\n' ...
        'Mean Overlap ± SD:\n' ...
        'ML: %.3f ± %.3f\n' ...
        'AP: %.3f ± %.3f\n' ...
        'DV: %.3f ± %.3f\n\n' ...
        'Strongest Effects:\n' ...
        'ML: %s (%.3f)\n' ...
        'AP: %s (%.3f)\n' ...
        'DV: %s (%.3f)\n\n' ...
        'Cell Types Compared:\n%s vs %s'], ...
        n_sig_ml, n_total, 100*n_sig_ml/n_total, ...
        n_sig_ap, n_total, 100*n_sig_ap/n_total, ...
        n_sig_dv, n_total, 100*n_sig_dv/n_total, ...
        mean(ML_overlap(valid_idx)), std(ML_overlap(valid_idx)), ...
        mean(AP_overlap(valid_idx)), std(AP_overlap(valid_idx)), ...
        mean(DV_overlap(valid_idx)), std(DV_overlap(valid_idx)), ...
        valid_areas{idx_ml}, max_ml, ...
        valid_areas{idx_ap}, max_ap, ...
        valid_areas{idx_dv}, max_dv, ...
        myType{1}, myType{2});
    
    text(0.05, 0.95, summary_text, 'VerticalAlignment', 'top', ...
         'FontSize', 10, 'Interpreter', 'tex');
end
%% preprocessing for topographical overlap
% ===== SETTINGS =====
use_weighted_mean = false;  % Set to false for simple (unweighted) mean
% ====================

% Initialize cell arrays to store results for each celltype
ML_ksd_all = {};
AP_ksd_all = {};
DV_ksd_all = {};
retinotopy_overlap_all = {};
numbers_all = {};

for c = 1:length(celltypes)
    ML_ksd = [];
    AP_ksd = [];
    DV_ksd = [];
    weights = [];  % Store weights for weighted mean
    
    % Get number of animals for this celltype
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    
    for a = 1:size(coord_fin, 2) % loop over areas
        for i = 1:num_animals
            n1 = size(coord_fin{c,a}{1,1}{i}, 1);
            n2 = size(coord_fin{c,a}{1,2}{i}, 1);
            
            %---ML---%
            if n1 >= 31 && n2 >= 31
                % Calculate weight (geometric mean of neuron counts)
                weights(a,i) = sqrt(n1 * n2);
                
                [y_med,~] = ksdensity(coord_fin{c,a}{1,1}{i}(:,3)); %ML
                x_med = linspace(min(coord_fin{c,a}{1,1}{i}(:,3)), max(coord_fin{c,a}{1,1}{i}(:,3)), 100);
                [y_lat,~] = ksdensity(coord_fin{c,a}{1,2}{i}(:,3));
                x_lat = linspace(min(coord_fin{c,a}{1,2}{i}(:,3)), max(coord_fin{c,a}{1,2}{i}(:,3)), 100);
                newx = linspace(min([min(coord_fin{c,a}{1,1}{i}(:,3)) min(coord_fin{c,a}{1,2}{i}(:,3))]), ...
                    max([max(coord_fin{c,a}{1,1}{i}(:,3)) max(coord_fin{c,a}{1,2}{i}(:,3))]), 100);
                new_y_med = interp1(x_med, y_med, newx);
                new_y_lat = interp1(x_lat, y_lat, newx);
                pmin = min(new_y_med, new_y_lat);
                pmax = max(new_y_med, new_y_lat);
                overlap = sum(pmin) ./ sum(pmax);
                ML_ksd(a,i) = overlap;
                
                %---AP---%
                [y_med,~] = ksdensity(coord_fin{c,a}{1,1}{i}(:,1)); %AP
                x_med = linspace(min(coord_fin{c,a}{1,1}{i}(:,1)), max(coord_fin{c,a}{1,1}{i}(:,1)), 100);
                [y_lat,~] = ksdensity(coord_fin{c,a}{1,2}{i}(:,1));
                x_lat = linspace(min(coord_fin{c,a}{1,2}{i}(:,1)), max(coord_fin{c,a}{1,2}{i}(:,1)), 100);
                newx = linspace(min([min(coord_fin{c,a}{1,1}{i}(:,1)) min(coord_fin{c,a}{1,2}{i}(:,1))]), ...
                    max([max(coord_fin{c,a}{1,1}{i}(:,1)) max(coord_fin{c,a}{1,2}{i}(:,1))]), 100);
                new_y_med = interp1(x_med, y_med, newx);
                new_y_lat = interp1(x_lat, y_lat, newx);
                pmin = min(new_y_med, new_y_lat);
                pmax = max(new_y_med, new_y_lat);
                overlap = sum(pmin) ./ sum(pmax);
                AP_ksd(a,i) = overlap;
                
                %---DV---%
                [y_med,~] = ksdensity(coord_fin{c,a}{1,1}{i}(:,2)); %DV
                x_med = linspace(min(coord_fin{c,a}{1,1}{i}(:,2)), max(coord_fin{c,a}{1,1}{i}(:,2)), 100);
                [y_lat,~] = ksdensity(coord_fin{c,a}{1,2}{i}(:,2));
                x_lat = linspace(min(coord_fin{c,a}{1,2}{i}(:,2)), max(coord_fin{c,a}{1,2}{i}(:,2)), 100);
                newx = linspace(min([min(coord_fin{c,a}{1,1}{i}(:,2)) min(coord_fin{c,a}{1,2}{i}(:,2))]), ...
                    max([max(coord_fin{c,a}{1,1}{i}(:,2)) max(coord_fin{c,a}{1,2}{i}(:,2))]), 100);
                new_y_med = interp1(x_med, y_med, newx);
                new_y_lat = interp1(x_lat, y_lat, newx);
                pmin = min(new_y_med, new_y_lat);
                pmax = max(new_y_med, new_y_lat);
                overlap = sum(pmin) ./ sum(pmax);
                DV_ksd(a,i) = overlap;
            else
                ML_ksd(a,i) = nan;
                AP_ksd(a,i) = nan;
                DV_ksd(a,i) = nan;
                weights(a,i) = nan;
            end
        end
    end
    
    % Apply mask
    mask = isnan(ML_ksd) | (ML_ksd == 0);
    ML_ksd(mask) = nan;
    AP_ksd(mask) = nan;
    DV_ksd(mask) = nan;
    weights(mask) = nan;
    
    % Compute means across animals (dimension 2)
    if use_weighted_mean
        % Weighted mean
        ML_ksd_mean = nan(size(ML_ksd, 1), 1);
        AP_ksd_mean = nan(size(AP_ksd, 1), 1);
        DV_ksd_mean = nan(size(DV_ksd, 1), 1);
        
        for a = 1:size(ML_ksd, 1)
            valid = ~isnan(ML_ksd(a,:));
            if any(valid)
                w = weights(a, valid);
                w_norm = w / sum(w);  % normalize weights to sum to 1
                ML_ksd_mean(a) = sum(ML_ksd(a, valid) .* w_norm);
                AP_ksd_mean(a) = sum(AP_ksd(a, valid) .* w_norm);
                DV_ksd_mean(a) = sum(DV_ksd(a, valid) .* w_norm);
            end
        end
    else
        % Simple unweighted mean
        ML_ksd_mean = mean(ML_ksd, 2, 'omitnan');
        AP_ksd_mean = mean(AP_ksd, 2, 'omitnan');
        DV_ksd_mean = mean(DV_ksd, 2, 'omitnan');
    end
    
    % Store results for this celltype
    ML_ksd_all{c} = ML_ksd_mean;
    AP_ksd_all{c} = AP_ksd_mean;
    DV_ksd_all{c} = DV_ksd_mean;
    
    % Create table for this celltype
    retinotopy_overlap_all{c} = table(areas, ML_ksd_mean, AP_ksd_mean, DV_ksd_mean, ...
        'VariableNames', {'areas', 'ML_ksd', 'AP_ksd', 'DV_ksd'});
    
    % Compute numbers (min across axes)
    numbers = [ML_ksd_mean'; AP_ksd_mean'; DV_ksd_mean'];
    numbers_all{c} = min(numbers, [], 1);
end

% Optional: Create a combined table with celltype as a variable
retinotopy_combined = table();
for c = 1:length(celltypes)
    temp = retinotopy_overlap_all{c};
    temp.celltype = repmat(celltypes(c), height(temp), 1);
    retinotopy_combined = [retinotopy_combined; temp];
end

% Print which method was used
if use_weighted_mean
    fprintf('Used WEIGHTED mean (by geometric mean of neuron counts)\n');
else
    fprintf('Used UNWEIGHTED (simple) mean\n');
end


%% Topographical Analysis - KDE version

% Topographical Analysis - KDE version
% Per-animal separation calculation with per-animal axis limits in plots

areas = string(union_labels(:,1));
n_areas = length(areas);
n_celltypes = length(celltypes);

% === CUSTOMIZE COLORS HERE ===
colors_ntsr = [col(3,:);   % upper
               col(4,:)];  % lower

colors_grp = [col(1,:);    % upper
              col(2,:)];   % lower

celltype_colors = {colors_ntsr, colors_grp};
% =============================

% Precompute grid of directions
n_theta = 18;
n_phi = 36;
theta_vals = linspace(0, pi, n_theta);
phi_vals = linspace(0, 2*pi, n_phi);
[THETA, PHI] = meshgrid(theta_vals, phi_vals);
THETA = THETA(:);
PHI = PHI(:);

DIR_X = sin(THETA) .* cos(PHI);
DIR_Y = sin(THETA) .* sin(PHI);
DIR_Z = cos(THETA);
all_directions = [DIR_X, DIR_Y, DIR_Z];
n_dir = size(all_directions, 1);

% Storage - per animal
separation_per_animal = cell(n_celltypes, n_areas);
direction1_per_animal = cell(n_celltypes, n_areas);
direction2_per_animal = cell(n_celltypes, n_areas);
sep2_per_animal = cell(n_celltypes, n_areas);

% Storage - mean across animals
max_separation_mean = cell(n_celltypes, n_areas);
second_separation_mean = cell(n_celltypes, n_areas);

for c = 1:n_celltypes
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    fprintf('\n=== Processing %s ===\n', celltypes{c});
    
    for a = 1:n_areas
        sep_animals = nan(num_animals, 1);
        sep2_animals = nan(num_animals, 1);
        dir1_animals = nan(num_animals, 3);
        dir2_animals = nan(num_animals, 3);
        
        for i = 1:num_animals
            upper_coords = coord_fin{c,a}{1,1}{i};
            lower_coords = coord_fin{c,a}{1,2}{i};
            
            n1 = size(upper_coords, 1);
            n2 = size(lower_coords, 1);
            
            if n1 < 31 || n2 < 31
                continue;
            end
            
            % Find optimal direction for THIS animal
            upper_proj = upper_coords * all_directions';
            lower_proj = lower_coords * all_directions';
            
            sep_vals = zeros(n_dir, 1);
            for d = 1:n_dir
                sep_vals(d) = calc_sep_kde(upper_proj(:,d), lower_proj(:,d));
            end
            
            [best_sep, best_idx] = max(sep_vals);
            best_dir = all_directions(best_idx, :)';
            
            sep_animals(i) = best_sep;
            dir1_animals(i,:) = best_dir';
            
            % Find second direction for THIS animal
            [second_dir, sep2] = find_second_dir_kde(upper_coords, lower_coords, best_dir);
            dir2_animals(i,:) = second_dir';
            sep2_animals(i) = sep2;
        end
        
        % Store per-animal results
        separation_per_animal{c,a} = sep_animals;
        direction1_per_animal{c,a} = dir1_animals;
        direction2_per_animal{c,a} = dir2_animals;
        sep2_per_animal{c,a} = sep2_animals;
        
        % Calculate means
        valid_idx = ~isnan(sep_animals);
        if sum(valid_idx) > 0
            max_separation_mean{c,a} = mean(sep_animals(valid_idx));
            second_separation_mean{c,a} = mean(sep2_animals(valid_idx));
            % max_separation_mean{c,a} = max(sep_animals(valid_idx));
            % second_separation_mean{c,a} = max(sep2_animals(valid_idx));
            
            fprintf('  %s: Sep1=%.3f, Sep2=%.3f (n=%d animals)\n', ...
                areas(a), max_separation_mean{c,a}, second_separation_mean{c,a}, sum(valid_idx));
        end
    end
end

fprintf('\nCalculations complete.\n');

%% Permutation null for separation metric
% Coarser grid for permutations only
% n_theta_null = 9;
% n_phi_null = 18;
n_theta_null = 18;
n_phi_null = 36;
theta_null = linspace(0, pi, n_theta_null);
phi_null = linspace(0, 2*pi, n_phi_null);
[TH_n, PH_n] = meshgrid(theta_null, phi_null);
TH_n = TH_n(:); PH_n = PH_n(:);
null_directions = [sin(TH_n).*cos(PH_n), sin(TH_n).*sin(PH_n), cos(TH_n)];
n_dir_null = size(null_directions, 1);


n_perm = 300;                      % 200 is enough if slow
null_mean_sep = nan(n_celltypes, n_areas, n_perm);
p_sep = nan(n_celltypes, n_areas);

for c = 1:n_celltypes
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    for a = 1:n_areas

        % --- collect animals that passed the 30-cell floor ---
        pooled = {}; n1_list = []; n2_list = [];
        for i = 1:num_animals
            u = coord_fin{c,a}{1,1}{i};
            l = coord_fin{c,a}{1,2}{i};
            if size(u,1) < 31 || size(l,1) < 31, continue; end
            pooled{end+1} = [u; l];          %#ok<SAGROW>
            n1_list(end+1) = size(u,1);      %#ok<SAGROW>
            n2_list(end+1) = size(l,1);      %#ok<SAGROW>
        end
        if isempty(pooled), continue; end
        n_upper_per_area{c,a} = n1_list;
        n_lower_per_area{c,a} = n2_list;

        % precompute projections once per animal (labels don't change these)
        proj_all = cellfun(@(X) X * null_directions', pooled, 'UniformOutput', false);

        nA = numel(pooled);
        perm_means = nan(n_perm, 1);

        parfor p = 1:n_perm                  % use 'for' if no toolbox
            sep_this_perm = nan(nA, 1);
            for k = 1:nA
                nTot = n1_list(k) + n2_list(k);
                idx  = randperm(nTot);
                iU   = idx(1:n1_list(k));
                iL   = idx(n1_list(k)+1:end);

                P = proj_all{k};
                sv = zeros(n_dir_null, 1);
                for d = 1:n_dir_null
                    sv(d) = calc_sep_kde(P(iU,d), P(iL,d));
                end
                sep_this_perm(k) = max(sv);
            end
            perm_means(p) = mean(sep_this_perm, 'omitnan');
        end

        null_mean_sep(c,a,:) = perm_means;

        obs = max_separation_mean{c,a};
        if ~isempty(obs) && ~isnan(obs)
            p_sep(c,a) = (1 + sum(perm_means >= obs)) / (n_perm + 1);
        end
    end
end

% FDR across areas, within celltype
p_sep_fdr = nan(size(p_sep));
for c = 1:n_celltypes
    v = p_sep(c,:); ok = ~isnan(v);
    p_sep_fdr(c,ok) = mafdr(v(ok), 'BHFDR', true);
end

null_p95 = prctile(null_mean_sep, 95, 3);          % n_celltypes x n_areas
obs_sep = max_separation_mean;
obs_sep(cellfun(@isempty, obs_sep)) = {NaN};
obs_sep = cell2mat(obs_sep);

excess = obs_sep - null_p95;

null_mean_val = mean(null_mean_sep, 3, 'omitnan');
null_std_val = std(null_mean_sep, 0, 3, 'omitnan');
effect_d = (obs_sep - null_mean_val) ./ null_std_val;

%% plot it (cohen's d)

effect_filtered = effect_d(:, ~both_nan);
effect_sorted = effect_filtered(:, sorted_indices);

cap = 40;
effect_capped = effect_sorted;
effect_capped(effect_capped > cap) = cap;

f2 = figure(2);
set(f2, 'Position', [73 200 992 180]);
hm2 = heatmap(effect_capped, 'Colormap', magma(300), ...
    'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 10);
hm2.Title = 'effect size (Cohen''s d relative to permutation null)';
hm2.MissingDataColor = [0.5 0.5 0.5];
hm2.XDisplayLabels = areas_sorted;
hm2.YDisplayLabels = celltypes;
hm2.ColorLimits = [0 cap];

s2 = struct(hm2);
s2.Axes.Title.FontWeight = 'normal';
s2.XAxis.TickLabelRotation = 90;
set(hm2, 'OuterPosition', [0 0 1 0.8]);
s2.Colorbar.Label.String = "Cohen's d";

c = 1; a = find(strcmp(areas, "APN"));
nd = squeeze(null_mean_sep(c,a,:));
obs = max_separation_mean{c,a};
null_mu = mean(nd, 'omitnan');
null_sd = std(nd, 'omitnan');
d = (obs - null_mu) / null_sd;

figure; histogram(nd, 30, 'EdgeColor', 'none'); hold on
xline(obs, 'r', 'LineWidth', 2);
xline(prctile(nd, 95), 'k--', 'LineWidth', 2);
xline(null_mu, 'b', 'LineWidth', 1.5);
xline(null_mu + null_sd, 'b--', 'LineWidth', 1);
xlim([0 1]);
legend('null', 'observed', 'null 95th', 'null mean', 'mean + 1 SD');
title(sprintf('APN — Cohen''s d = %.1f', d));

%% plot it
% excess=obs_sep - null_p95;
excess = obs_sep - null_mean_val;

excess_filtered = excess(:, ~both_nan);
excess_sorted = excess_filtered(:, sorted_indices);

f2 = figure(2);
set(f2, 'Position', [73 200 992 180]);

excess_clipped = excess_sorted;
excess_clipped(excess_clipped < 0) = 0;

lim = max(abs(excess_sorted(:)), [], 'omitnan');
hm2 = heatmap(excess_clipped, 'Colormap', magma(300), ...
    'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 10);
hm2.MissingDataColor = [0.5 0.5 0.5];
hm2.XDisplayLabels = areas_sorted;
hm2.YDisplayLabels = celltypes;
hm2.ColorLimits = [0 max(excess_clipped(:), [], 'omitnan')];
% hm2.ColorLimits = [-lim lim];



s2 = struct(hm2);
s2.Axes.Title.FontWeight = 'normal';
s2.XAxis.TickLabelRotation = 90;
set(hm2, 'OuterPosition', [0 0 1 0.8]);
s2.Colorbar.Label.String = ["excess"; "separation"];

c = 1; a = find(strcmp(areas, "VISp"));
figure; histogram(squeeze(null_mean_sep(c,a,:)), 30, 'Edgecolor', 'none'); hold on
xline(max_separation_mean{c,a}, 'r', 'LineWidth', 2);
% xline(prctile(squeeze(null_mean_sep(c,a,:)), 95), 'k--', 'LineWidth', 2); %for plotting the 95th percentile
xline(mean(squeeze(null_mean_sep(c,a,:)), 'omitnan'), 'b', 'LineWidth', 2); %for plotting the mean
legend('null','observed','null 95th');
xlim([0 1]);

rows = {};
for c = 1:n_celltypes
    for a = 1:n_areas
        if ~isempty(n_upper_per_area{c,a})
            rows{end+1} = {celltypes{c}, areas{a}, ...
                numel(n_upper_per_area{c,a}), ...
                median(n_upper_per_area{c,a}), ...
                median(n_lower_per_area{c,a}), ...
                sum(n_upper_per_area{c,a}), ...
                sum(n_lower_per_area{c,a})};
        end
    end
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'CellType','Area','nAnimals', ...
     'MedianUpper','MedianLower', ...
     'TotalUpper','TotalLower'});
disp(T);
writetable(T, 'separation_neuron_counts.xlsx');

%% check values
c = 1;
for a = 1:n_areas
    nd = squeeze(null_mean_sep(c,a,:));
    if all(isnan(nd)) || isempty(max_separation_mean{c,a}), continue; end
    fprintf('%-8s obs=%.3f  null_med=%.3f  null95=%.3f  z=%+.2f\n', ...
        areas(a), max_separation_mean{c,a}, median(nd), prctile(nd,95), ...
        (max_separation_mean{c,a} - median(nd)) / std(nd));
end

%% Check axis consistency across animals
axis_spread = nan(n_celltypes, n_areas);

for c = 1:n_celltypes
    for a = 1:n_areas
        dirs = direction1_per_animal{c,a};
        valid = ~any(isnan(dirs), 2);
        dirs = dirs(valid, :);
        n = size(dirs, 1);
        if n < 3, continue; end
        
        % align: flip vectors pointing away from the mean
        mean_dir = mean(dirs, 1);
        mean_dir = mean_dir / norm(mean_dir);
        for i = 1:n
            if dot(dirs(i,:), mean_dir) < 0
                dirs(i,:) = -dirs(i,:);
            end
        end
        
        % recompute mean after alignment
        mean_dir = mean(dirs, 1);
        mean_dir = mean_dir / norm(mean_dir);
        
        % angular deviation of each animal from the mean
        angles = acosd(min(1, max(-1, dirs * mean_dir')));
        axis_spread(c,a) = median(angles);
        
        fprintf('%s %-8s: median deviation = %.1f deg (n=%d animals)\n', ...
            celltypes{c}, areas(a), median(angles), n);
    end
end

%% Export axis consistency table
rows = {};
for c = 1:n_celltypes
    for a = 1:n_areas
        dirs = direction1_per_animal{c,a};
        valid = ~any(isnan(dirs), 2);
        dirs = dirs(valid, :);
        n = size(dirs, 1);
        if n < 1, continue; end
        
        % axis deviation
        if n >= 2
            mean_dir = mean(dirs, 1);
            mean_dir = mean_dir / norm(mean_dir);
            for i = 1:n
                if dot(dirs(i,:), mean_dir) < 0
                    dirs(i,:) = -dirs(i,:);
                end
            end
            mean_dir = mean(dirs, 1);
            mean_dir = mean_dir / norm(mean_dir);
            angles = acosd(min(1, max(-1, dirs * mean_dir')));
            med_dev = median(angles);
        else
            med_dev = NaN;
        end
        
        % observed separation
        if ~isempty(max_separation_mean{c,a})
            obs = max_separation_mean{c,a};
        else
            obs = NaN;
        end
        
        % null values
        nd = squeeze(null_mean_sep(c,a,:));
        if all(isnan(nd))
            null_med = NaN;
            null_95 = NaN;
        else
            null_med = median(nd, 'omitnan');
            null_95 = prctile(nd, 95);
        end
        
        % neuron counts
        if ~isempty(n_upper_per_area{c,a})
            n_up = sum(n_upper_per_area{c,a});
            n_low = sum(n_lower_per_area{c,a});
        else
            n_up = NaN;
            n_low = NaN;
        end
        
        rows{end+1} = {celltypes{c}, areas{a}, n, n_up, n_low, ...
            round(obs,3), round(null_med,3), round(null_95,3), ...
            round(obs - null_95, 3), round(med_dev,1)};
    end
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'CellType','Area','nAnimals','nUpper','nLower', ...
     'ObsSeparation','NullMedian','Null95th','Excess','MedianAxisDev_deg'});

disp(T);
writetable(T, 'separation_metric_table.xlsx');

%% Check how peaked the separation-vs-direction curve is
for c = 1:n_celltypes
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    for a = 1:n_areas
        ratios = [];
        for i = 1:num_animals
            upper_coords = coord_fin{c,a}{1,1}{i};
            lower_coords = coord_fin{c,a}{1,2}{i};
            if size(upper_coords,1) < 31 || size(lower_coords,1) < 31, continue; end
            
            upper_proj = upper_coords * all_directions';
            lower_proj = lower_coords * all_directions';
            sv = zeros(n_dir, 1);
            for d = 1:n_dir
                sv(d) = calc_sep_kde(upper_proj(:,d), lower_proj(:,d));
            end
            ratios(end+1) = max(sv) / mean(sv);
        end
        if ~isempty(ratios)
            fprintf('%s %-8s: max/mean = %.2f (n=%d)\n', ...
                celltypes{c}, areas(a), mean(ratios), numel(ratios));
        end
    end
end
%% ========================================================================
% POOLED 2D SCATTER PLOTS - ALL AREAS
% Axis 1: optimal separation direction
% Axis 2: principal component of residual variance (perpendicular to axis 1)
% =========================================================================

for c = 1:n_celltypes
    valid_areas = [];
    for a = 1:n_areas
        if ~isempty(max_separation_mean{c,a}) && ~isnan(max_separation_mean{c,a})
            valid_areas = [valid_areas, a];
        end
    end
    if isempty(valid_areas)
        continue;
    end
    
    colors = celltype_colors{c};
    color_upper = colors(1,:);
    color_lower = colors(2,:);
    
    figure('Position', [100, 100, 1600, 1200]);
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('%s neurons - Optimal separation (pooled)', celltypes{c}), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    
    for idx = 1:length(valid_areas)
        a = valid_areas(idx);
        nexttile;
        hold on;
        
        % Pool data across animals
        upper_all = [];
        lower_all = [];
        for i = 1:num_animals
            if ~isempty(coord_fin{c,a}{1,1}{i})
                upper_all = [upper_all; coord_fin{c,a}{1,1}{i}];
            end
            if ~isempty(coord_fin{c,a}{1,2}{i})
                lower_all = [lower_all; coord_fin{c,a}{1,2}{i}];
            end
        end
        
        if isempty(upper_all) || isempty(lower_all)
            continue;
        end
        
        % Compute optimal direction on POOLED data
        upper_proj_all = upper_all * all_directions';
        lower_proj_all = lower_all * all_directions';
        sep_vals = zeros(n_dir, 1);
        for d = 1:n_dir
            sep_vals(d) = calc_sep_kde(upper_proj_all(:,d), lower_proj_all(:,d));
        end
        [~, best_idx] = max(sep_vals);
        dir1 = all_directions(best_idx, :)';
        dir1 = dir1 / norm(dir1);
        
        % Axis 2: PCA on the plane perpendicular to dir1
        all_coords = [upper_all; lower_all];
        % Project out dir1 component
        proj_on_dir1 = all_coords * dir1;  % Nx1
        residuals = all_coords - proj_on_dir1 * dir1';  % Nx3, orthogonal to dir1
        
        % PCA on residuals to find direction of max spread
        [V, ~] = eig(cov(residuals));
        % Eigenvectors sorted ascending by eigenvalue, take last (largest)
        dir2 = V(:, end);
        dir2 = dir2 / norm(dir2);
        
        % Ensure dir2 is truly orthogonal (numerical cleanup)
        dir2 = dir2 - dot(dir2, dir1) * dir1;
        dir2 = dir2 / norm(dir2);
        
        % Project pooled data
        upper_x = upper_all * dir1;
        upper_y = upper_all * dir2;
        lower_x = lower_all * dir1;
        lower_y = lower_all * dir2;
        
        % Scatter plot
        s_upper = scatter(upper_x, upper_y, 8, color_upper, 'filled', 'MarkerFaceAlpha', 0.45);
        s_lower = scatter(lower_x, lower_y, 8, color_lower, 'filled', 'MarkerFaceAlpha', 0.45);
        uistack(s_lower, 'bottom');
        
        % Centroids
        scatter(mean(upper_x), mean(upper_y), 100, color_upper, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        scatter(mean(lower_x), mean(lower_y), 100, color_lower, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        xlabel('Separation axis');
        ylabel('Max spread (perp.)');
        title(sprintf('%s (sep=%.2f)', areas(a), max_separation_mean{c,a}), 'Interpreter', 'none');
        axis equal;
        grid on;
        box on;
        
        if idx == 1
            legend({'Upper', 'Lower'}, 'Location', 'best');
        end
    end
    set(gcf, 'Renderer', 'painters');
end

%% ========================================================================
% MEAN 1D DENSITY PLOTS - ALL AREAS
% Shows mean density across animals along optimal axis
% =========================================================================

for c = 1:n_celltypes
    valid_areas = [];
    for a = 1:n_areas
        if ~isempty(max_separation_mean{c,a}) && ~isnan(max_separation_mean{c,a})
            valid_areas = [valid_areas, a];
        end
    end
    
    if isempty(valid_areas)
        continue;
    end
    
    colors = celltype_colors{c};
    color_upper = colors(1,:);
    color_lower = colors(2,:);
    
    figure('Position', [100, 100, 1600, 1200]);
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('%s neurons - Mean 1D density along optimal axis', celltypes{c}), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    
    for idx = 1:length(valid_areas)
        a = valid_areas(idx);
        
        nexttile;
        hold on;
        
        % Get valid animals
        valid_idx = find(~isnan(separation_per_animal{c,a}));
        
        if isempty(valid_idx)
            continue;
        end
        
        % Collect all projected data to determine common x range
        all_proj = [];
        for i = valid_idx'
            upper_coords = coord_fin{c,a}{1,1}{i};
            lower_coords = coord_fin{c,a}{1,2}{i};
            dir1 = direction1_per_animal{c,a}(i,:)';
            dir1 = dir1 / norm(dir1);
            all_proj = [all_proj; upper_coords * dir1; lower_coords * dir1];
        end
        
        x_range = linspace(min(all_proj)-0.1, max(all_proj)+0.1, 300);
        
        % Compute KDE for each animal and average
        kde_upper_all = zeros(length(valid_idx), length(x_range));
        kde_lower_all = zeros(length(valid_idx), length(x_range));
        
        for ii = 1:length(valid_idx)
            i = valid_idx(ii);
            upper_coords = coord_fin{c,a}{1,1}{i};
            lower_coords = coord_fin{c,a}{1,2}{i};
            
            dir1 = direction1_per_animal{c,a}(i,:)';
            dir1 = dir1 / norm(dir1);
            
            upper_proj = upper_coords * dir1;
            lower_proj = lower_coords * dir1;
            
            % KDE with support
            support_upper = [min(upper_proj)-0.1, max(upper_proj)+0.1];
            support_lower = [min(lower_proj)-0.1, max(lower_proj)+0.1];
            
            [kde_u, x_u] = ksdensity(upper_proj, 'NumPoints', 300, 'Support', support_upper);
            [kde_l, x_l] = ksdensity(lower_proj, 'NumPoints', 300, 'Support', support_lower);
            
            % Interpolate to common x range
            kde_upper_all(ii,:) = interp1(x_u, kde_u, x_range, 'linear', 0);
            kde_lower_all(ii,:) = interp1(x_l, kde_l, x_range, 'linear', 0);
        end
        
        % Mean across animals
        kde_upper_mean = mean(kde_upper_all, 1);
        kde_lower_mean = mean(kde_lower_all, 1);
        
        % Plot
        plot(x_range, kde_upper_mean, '-', 'Color', color_upper, 'LineWidth', 2);
        plot(x_range, kde_lower_mean, '-', 'Color', color_lower, 'LineWidth', 2);
        
        % Shade overlap
        overlap_density = min(kde_upper_mean, kde_lower_mean);
        fill([x_range, fliplr(x_range)], ...
             [overlap_density, zeros(size(overlap_density))], ...
             [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        xlabel('Position along optimal axis (mm)');
        ylabel('Density');
        title(sprintf('%s (sep=%.2f)', areas(a), max_separation_mean{c,a}), 'Interpreter', 'none');
        grid on;
        box on;
        
        if idx == 1
            legend({'Upper', 'Lower', 'Overlap'}, 'Location', 'best');
        end
    end
end


%% ========================================================================
% CUSTOM AREA PLOTTING SECTION
% Change the area name below and run this section to plot individual animals
% =========================================================================

plot_area = 'PBG';  % <-- CHANGE THIS TO THE AREA YOU WANT TO PLOT

% Find area index
area_idx = find(strcmpi(areas, plot_area));
if isempty(area_idx)
    error('Area "%s" not found. Available areas: %s', plot_area, strjoin(areas, ', '));
end

for c = 1:n_celltypes
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    
    colors = celltype_colors{c};
    color_upper = colors(1,:);
    color_lower = colors(2,:);
    
    % Count valid animals for this area
    sep_vals = separation_per_animal{c, area_idx};
    valid_animals = find(~isnan(sep_vals));
    n_valid = length(valid_animals);
    
    if n_valid == 0
        fprintf('No valid animals for %s in %s\n', celltypes{c}, plot_area);
        continue;
    end
    
    % Create figure
    f = figure();
    bigT = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(bigT, sprintf('%s - %s', celltypes{c}, plot_area), 'FontSize', 14, 'FontWeight', 'bold');
    
    counter = 1;
    for fi = valid_animals'
        upper_coords = coord_fin{c, area_idx}{1,1}{fi};
        lower_coords = coord_fin{c, area_idx}{1,2}{fi};
        
        % Get this animal's optimal directions
        dir1 = direction1_per_animal{c, area_idx}(fi,:)';
        dir2 = direction2_per_animal{c, area_idx}(fi,:)';
        
        dir1 = dir1 / norm(dir1);
        dir2 = dir2 - dot(dir2, dir1) * dir1;
        dir2 = dir2 / norm(dir2);
        
        % Project onto optimal axes
        upper_x = upper_coords * dir1;
        upper_y = upper_coords * dir2;
        lower_x = lower_coords * dir1;
        lower_y = lower_coords * dir2;
        
        % Per-animal axis limits with EQUAL SCALE
        all_x = [upper_x; lower_x];
        all_y = [upper_y; lower_y];
        x_range = max(all_x) - min(all_x);
        y_range = max(all_y) - min(all_y);
        max_range = max(x_range, y_range);
        
        x_center = (max(all_x) + min(all_x)) / 2;
        y_center = (max(all_y) + min(all_y)) / 2;
        
        ax1_lim = [x_center - max_range/2 - 0.1, x_center + max_range/2 + 0.1];
        ax2_lim = [y_center - max_range/2 - 0.1, y_center + max_range/2 + 0.1];
        
        sep1 = separation_per_animal{c, area_idx}(fi);
        sep2 = sep2_per_animal{c, area_idx}(fi);
        
        % Create nested tiledlayout
        nexttile(bigT);
        set(gca, 'box', 'off', 'visible', 'off');
        
        fact = 2;
        til = tiledlayout(bigT, 3*fact, 3*fact, 'TileSpacing', 'none', 'Padding', 'tight');
        til.Layout.Tile = counter;
        til.Title.String = sprintf('Animal %d (sep=%.2f)', fi, sep1);
        til.Title.FontSize = 12;
        til.Title.FontWeight = 'bold';
        
        % Horizontal density curve (top)
        nexttile(til, [1, 2*fact+1]);
        support_upper_x = [min(upper_x)-0.1, max(upper_x)+0.1];
        support_lower_x = [min(lower_x)-0.1, max(lower_x)+0.1];
        [y_upper, x_upper] = ksdensity(upper_x, 'NumPoints', 300, 'Support', support_upper_x);
        plot(x_upper, y_upper, 'LineWidth', 2, 'Color', color_upper); hold on;
        [y_lower, x_lower] = ksdensity(lower_x, 'NumPoints', 300, 'Support', support_lower_x);
        plot(x_lower, y_lower, 'LineWidth', 2, 'Color', color_lower);
        xlim(gca, ax1_lim);
        set(gca, 'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
        
        % Scatter plot (center)
        nexttile(til, 4*fact-1, [2*fact+1, 2*fact+1]);
        s_upper = scatter(upper_x, upper_y, 10, 'filled', 'CData', color_upper, 'MarkerFaceAlpha', 0.6); hold on;
        s_lower = scatter(lower_x, lower_y, 10, 'filled', 'CData', color_lower, 'MarkerFaceAlpha', 0.6);
        uistack(s_lower, 'bottom');
        
        set(gca, 'box', 'on', 'Color', 'none');
        xlim(gca, ax1_lim);
        ylim(gca, ax2_lim);
        % axis equal;
        xlim(gca, ax1_lim);
        ylim(gca, ax2_lim);
        
        xlabel('Axis 1', 'FontSize', 12);
        ylabel('Axis 2', 'FontSize', 12);
        set(gca, 'FontSize', 10);
        
        % Vertical density curve (right)
        nexttile(til, 6*fact, [2*fact+1, 1]);
        support_upper_y = [min(upper_y)-0.1, max(upper_y)+0.1];
        support_lower_y = [min(lower_y)-0.1, max(lower_y)+0.1];
        [y_upper_v, x_upper_v] = ksdensity(upper_y, 'NumPoints', 300, 'Support', support_upper_y);
        plot(y_upper_v, x_upper_v, 'LineWidth', 2, 'Color', color_upper); hold on;
        [y_lower_v, x_lower_v] = ksdensity(lower_y, 'NumPoints', 300, 'Support', support_lower_y);
        plot(y_lower_v, x_lower_v, 'LineWidth', 2, 'Color', color_lower);
        ylim(gca, ax2_lim);
        set(gca, 'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
        
        counter = counter + 1;
    end
    
    % Set figure size
    fig_width = 250 * n_valid;
    set(f, 'Position', [29, 150, fig_width, fig_width * 0.8], 'Renderer', 'Painters');
end


% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function sep = calc_sep_kde(x1, x2)
    all_x = [x1; x2];
    x_range = linspace(min(all_x), max(all_x), 100);
    
    if length(x1) < 2 || length(x2) < 2
        sep = NaN;
        return;
    end
    
    kde1 = ksdensity(x1, x_range);
    kde2 = ksdensity(x2, x_range);
    
    dx = x_range(2) - x_range(1);
    kde1 = kde1 / (sum(kde1) * dx);
    kde2 = kde2 / (sum(kde2) * dx);
    
    overlap = sum(min(kde1, kde2)) * dx;
    sep = 1 - overlap;
end

function [second_dir, second_sep] = find_second_dir_kde(upper_coords, lower_coords, first_dir)
    first_dir = first_dir / norm(first_dir);
    
    if abs(first_dir(1)) < 0.9
        v1 = cross(first_dir, [1;0;0]);
    else
        v1 = cross(first_dir, [0;1;0]);
    end
    v1 = v1 / norm(v1);
    v2 = cross(first_dir, v1);
    v2 = v2 / norm(v2);
    
    angles = linspace(0, pi, 19);
    best_sep = 0;
    best_angle = 0;
    
    for angle = angles
        test_dir = cos(angle) * v1 + sin(angle) * v2;
        upper_proj = upper_coords * test_dir;
        lower_proj = lower_coords * test_dir;
        sep = calc_sep_kde(upper_proj, lower_proj);
        
        if sep > best_sep
            best_sep = sep;
            best_angle = angle;
        end
    end
    
    second_dir = cos(best_angle) * v1 + sin(best_angle) * v2;
    second_sep = best_sep;
end

%% Topographical Analysis - LDA version
% Uses Fisher's Linear Discriminant Analysis
% Separation metric: Fisher criterion (between-class / within-class variance)

areas = string(union_labels(:,1));
n_areas = length(areas);
n_celltypes = length(celltypes);

% === CUSTOMIZE COLORS HERE ===
colors_ntsr = [col(3,:);   % upper
               col(4,:)];  % lower

colors_grp = [col(1,:);    % upper
              col(2,:)];   % lower

celltype_colors = {colors_ntsr, colors_grp};
% =============================

% Storage - per animal
separation_per_animal = cell(n_celltypes, n_areas);
direction1_per_animal = cell(n_celltypes, n_areas);
direction2_per_animal = cell(n_celltypes, n_areas);
sep2_per_animal = cell(n_celltypes, n_areas);

% Storage - mean across animals
max_separation_mean = cell(n_celltypes, n_areas);
second_separation_mean = cell(n_celltypes, n_areas);

for c = 1:n_celltypes
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    fprintf('\n=== Processing %s (LDA) ===\n', celltypes{c});
    
    for a = 1:n_areas
        sep_animals = nan(num_animals, 1);
        sep2_animals = nan(num_animals, 1);
        dir1_animals = nan(num_animals, 3);
        dir2_animals = nan(num_animals, 3);
        
        for i = 1:num_animals
            upper_coords = coord_fin{c,a}{1,1}{i};
            lower_coords = coord_fin{c,a}{1,2}{i};
            
            n1 = size(upper_coords, 1);
            n2 = size(lower_coords, 1);
            
            if n1 < 30 || n2 < 30
                continue;
            end
            
            % Fisher's LDA
            mean_upper = mean(upper_coords, 1)';
            mean_lower = mean(lower_coords, 1)';
            
            % Within-class scatter
            upper_centered = upper_coords - mean_upper';
            lower_centered = lower_coords - mean_lower';
            S_within = (upper_centered' * upper_centered + lower_centered' * lower_centered) / (n1 + n2 - 2);
            
            % Regularization for stability
            reg = 1e-4 * trace(S_within) / 3;
            S_within_reg = S_within + reg * eye(3);
            
            % LDA direction: S_within^(-1) * (mean_upper - mean_lower)
            mean_diff = mean_upper - mean_lower;
            lda_dir = S_within_reg \ mean_diff;
            lda_dir = lda_dir / norm(lda_dir);
            
            % Fisher criterion as separation metric
            % d' (d-prime): distance between means / pooled std
            upper_proj = upper_coords * lda_dir;
            lower_proj = lower_coords * lda_dir;
            d_prime = abs(mean(upper_proj) - mean(lower_proj)) / sqrt(0.5 * (var(upper_proj) + var(lower_proj)));
            
            sep_animals(i) = d_prime;
            dir1_animals(i,:) = lda_dir';
            
            % Find second direction (orthogonal to LDA direction, maximize Fisher criterion)
            [second_dir, sep2] = find_second_dir_lda(upper_coords, lower_coords, lda_dir);
            dir2_animals(i,:) = second_dir';
            sep2_animals(i) = sep2;
        end
        
        % Store per-animal results
        separation_per_animal{c,a} = sep_animals;
        direction1_per_animal{c,a} = dir1_animals;
        direction2_per_animal{c,a} = dir2_animals;
        sep2_per_animal{c,a} = sep2_animals;
        
        % Calculate means
        valid_idx = ~isnan(sep_animals);
        if sum(valid_idx) > 0
            max_separation_mean{c,a} = mean(sep_animals(valid_idx));
            second_separation_mean{c,a} = mean(sep2_animals(valid_idx));
            
            fprintf('  %s: d''1=%.3f, d''2=%.3f (n=%d animals)\n', ...
                areas(a), max_separation_mean{c,a}, second_separation_mean{c,a}, sum(valid_idx));
        end
    end
end

fprintf('\nLDA calculations complete.\n');

%% ========================================================================
% POOLED 2D SCATTER PLOTS - ALL AREAS - LDA
% =========================================================================

for c = 1:n_celltypes
    valid_areas = [];
    for a = 1:n_areas
        if ~isempty(max_separation_mean{c,a}) && ~isnan(max_separation_mean{c,a})
            valid_areas = [valid_areas, a];
        end
    end
    
    if isempty(valid_areas)
        continue;
    end
    
    colors = celltype_colors{c};
    color_upper = colors(1,:);
    color_lower = colors(2,:);
    
    figure('Position', [100, 100, 1600, 1200]);
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('%s neurons - LDA separation (pooled)', celltypes{c}), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    
    for idx = 1:length(valid_areas)
        a = valid_areas(idx);
        
        nexttile;
        hold on;
        
        % Pool data across animals
        upper_all = [];
        lower_all = [];
        for i = 1:num_animals
            if ~isempty(coord_fin{c,a}{1,1}{i})
                upper_all = [upper_all; coord_fin{c,a}{1,1}{i}];
            end
            if ~isempty(coord_fin{c,a}{1,2}{i})
                lower_all = [lower_all; coord_fin{c,a}{1,2}{i}];
            end
        end
        
        if isempty(upper_all) || isempty(lower_all)
            continue;
        end
        
        % Compute LDA directions on POOLED data
        n1 = size(upper_all, 1);
        n2 = size(lower_all, 1);
        
        mean_upper = mean(upper_all, 1)';
        mean_lower = mean(lower_all, 1)';
        
        upper_centered = upper_all - mean_upper';
        lower_centered = lower_all - mean_lower';
        S_within = (upper_centered' * upper_centered + lower_centered' * lower_centered) / (n1 + n2 - 2);
        
        reg = 1e-4 * trace(S_within) / 3;
        S_within_reg = S_within + reg * eye(3);
        
        mean_diff = mean_upper - mean_lower;
        dir1 = S_within_reg \ mean_diff;
        dir1 = dir1 / norm(dir1);
        
        % Find second LDA direction on pooled data
        [dir2, ~] = find_second_dir_lda(upper_all, lower_all, dir1);
        
        dir2 = dir2 - dot(dir2, dir1) * dir1;
        dir2 = dir2 / norm(dir2);
        
        % Project pooled data
        upper_x = upper_all * dir1;
        upper_y = upper_all * dir2;
        lower_x = lower_all * dir1;
        lower_y = lower_all * dir2;
        
        % Scatter plot
        s_upper = scatter(upper_x, upper_y, 8, color_upper, 'filled', 'MarkerFaceAlpha', 0.3);
        s_lower = scatter(lower_x, lower_y, 8, color_lower, 'filled', 'MarkerFaceAlpha', 0.3);
        uistack(s_lower, 'bottom');
        
        % Centroids
        scatter(mean(upper_x), mean(upper_y), 100, color_upper, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        scatter(mean(lower_x), mean(lower_y), 100, color_lower, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        xlabel('LDA 1');
        ylabel('LDA 2');
        title(sprintf('%s (d''=%.2f)', areas(a), max_separation_mean{c,a}), 'Interpreter', 'none');
        axis equal;
        grid on;
        box on;
        
        if idx == 1
            legend({'Upper', 'Lower'}, 'Location', 'best');
        end
    end
end


%% ========================================================================
% CUSTOM AREA PLOTTING SECTION
% Change the area name below and run this section to plot individual animals
% =========================================================================

plot_area = 'CUN';  % <-- CHANGE THIS TO THE AREA YOU WANT TO PLOT

% Find area index
area_idx = find(strcmpi(areas, plot_area));
if isempty(area_idx)
    error('Area "%s" not found. Available areas: %s', plot_area, strjoin(areas, ', '));
end

for c = 1:n_celltypes
    num_animals = size(coord_fin{c,1}{1,1}, 2);
    
    colors = celltype_colors{c};
    color_upper = colors(1,:);
    color_lower = colors(2,:);
    
    % Count valid animals for this area
    sep_vals = separation_per_animal{c, area_idx};
    valid_animals = find(~isnan(sep_vals));
    n_valid = length(valid_animals);
    
    if n_valid == 0
        fprintf('No valid animals for %s in %s\n', celltypes{c}, plot_area);
        continue;
    end
    
    % Create figure
    f = figure();
    bigT = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(bigT, sprintf('%s - %s (LDA)', celltypes{c}, plot_area), 'FontSize', 14, 'FontWeight', 'bold');
    
    counter = 1;
    for fi = valid_animals'
        upper_coords = coord_fin{c, area_idx}{1,1}{fi};
        lower_coords = coord_fin{c, area_idx}{1,2}{fi};
        
        % Get this animal's LDA directions
        dir1 = direction1_per_animal{c, area_idx}(fi,:)';
        dir2 = direction2_per_animal{c, area_idx}(fi,:)';
        
        dir1 = dir1 / norm(dir1);
        dir2 = dir2 - dot(dir2, dir1) * dir1;
        dir2 = dir2 / norm(dir2);
        
        % Project onto LDA axes
        upper_x = upper_coords * dir1;
        upper_y = upper_coords * dir2;
        lower_x = lower_coords * dir1;
        lower_y = lower_coords * dir2;
        
        % Per-animal axis limits with EQUAL SCALE
        all_x = [upper_x; lower_x];
        all_y = [upper_y; lower_y];
        x_range = max(all_x) - min(all_x);
        y_range = max(all_y) - min(all_y);
        max_range = max(x_range, y_range);
        
        x_center = (max(all_x) + min(all_x)) / 2;
        y_center = (max(all_y) + min(all_y)) / 2;
        
        ax1_lim = [x_center - max_range/2 - 0.1, x_center + max_range/2 + 0.1];
        ax2_lim = [y_center - max_range/2 - 0.1, y_center + max_range/2 + 0.1];
        
        sep1 = separation_per_animal{c, area_idx}(fi);
        sep2 = sep2_per_animal{c, area_idx}(fi);
        
        % Create nested tiledlayout
        nexttile(bigT);
        set(gca, 'box', 'off', 'visible', 'off');
        
        fact = 2;
        til = tiledlayout(bigT, 3*fact, 3*fact, 'TileSpacing', 'none', 'Padding', 'tight');
        til.Layout.Tile = counter;
        til.Title.String = sprintf('Animal %d (d''=%.2f)', fi, sep1);
        til.Title.FontSize = 12;
        til.Title.FontWeight = 'bold';
        
        % Horizontal density curve (top)
        nexttile(til, [1, 2*fact+1]);
        support_upper_x = [min(upper_x)-0.1, max(upper_x)+0.1];
        support_lower_x = [min(lower_x)-0.1, max(lower_x)+0.1];
        [y_upper, x_upper] = ksdensity(upper_x, 'NumPoints', 300, 'Support', support_upper_x);
        plot(x_upper, y_upper, 'LineWidth', 2, 'Color', color_upper); hold on;
        [y_lower, x_lower] = ksdensity(lower_x, 'NumPoints', 300, 'Support', support_lower_x);
        plot(x_lower, y_lower, 'LineWidth', 2, 'Color', color_lower);
        xlim(gca, ax1_lim);
        set(gca, 'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
        
        % Scatter plot (center)
        nexttile(til, 4*fact-1, [2*fact+1, 2*fact+1]);
        s_upper = scatter(upper_x, upper_y, 10, 'filled', 'CData', color_upper, 'MarkerFaceAlpha', 0.6); hold on;
        s_lower = scatter(lower_x, lower_y, 10, 'filled', 'CData', color_lower, 'MarkerFaceAlpha', 0.6);
        uistack(s_lower, 'bottom');
        
        set(gca, 'box', 'on', 'Color', 'none');
        xlim(gca, ax1_lim);
        ylim(gca, ax2_lim);
        % axis equal;
        xlim(gca, ax1_lim);
        ylim(gca, ax2_lim);
        
        xlabel('LDA 1', 'FontSize', 12);
        ylabel('LDA 2', 'FontSize', 12);
        set(gca, 'FontSize', 10);
        
        % Vertical density curve (right)
        nexttile(til, 6*fact, [2*fact+1, 1]);
        support_upper_y = [min(upper_y)-0.1, max(upper_y)+0.1];
        support_lower_y = [min(lower_y)-0.1, max(lower_y)+0.1];
        [y_upper_v, x_upper_v] = ksdensity(upper_y, 'NumPoints', 300, 'Support', support_upper_y);
        plot(y_upper_v, x_upper_v, 'LineWidth', 2, 'Color', color_upper); hold on;
        [y_lower_v, x_lower_v] = ksdensity(lower_y, 'NumPoints', 300, 'Support', support_lower_y);
        plot(y_lower_v, x_lower_v, 'LineWidth', 2, 'Color', color_lower);
        ylim(gca, ax2_lim);
        set(gca, 'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
        
        counter = counter + 1;
    end
    
    % Set figure size
    fig_width = 250 * n_valid;
    set(f, 'Position', [29, 150, fig_width, fig_width * 0.8], 'Renderer', 'Painters');
end


% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function [second_dir, second_dprime] = find_second_dir_lda(upper_coords, lower_coords, first_dir)
    % Find second LDA direction in the plane orthogonal to first_dir
    first_dir = first_dir / norm(first_dir);
    
    % Create orthonormal basis for the orthogonal plane
    if abs(first_dir(1)) < 0.9
        v1 = cross(first_dir, [1;0;0]);
    else
        v1 = cross(first_dir, [0;1;0]);
    end
    v1 = v1 / norm(v1);
    v2 = cross(first_dir, v1);
    v2 = v2 / norm(v2);
    
    % Project data onto orthogonal plane
    upper_2d = [upper_coords * v1, upper_coords * v2];
    lower_2d = [lower_coords * v1, lower_coords * v2];
    
    n1 = size(upper_2d, 1);
    n2 = size(lower_2d, 1);
    
    % 2D LDA in the orthogonal plane
    mean_upper_2d = mean(upper_2d, 1)';
    mean_lower_2d = mean(lower_2d, 1)';
    
    upper_centered = upper_2d - mean_upper_2d';
    lower_centered = lower_2d - mean_lower_2d';
    S_within_2d = (upper_centered' * upper_centered + lower_centered' * lower_centered) / (n1 + n2 - 2);
    
    reg = 1e-4 * trace(S_within_2d) / 2;
    S_within_2d_reg = S_within_2d + reg * eye(2);
    
    mean_diff_2d = mean_upper_2d - mean_lower_2d;
    lda_dir_2d = S_within_2d_reg \ mean_diff_2d;
    lda_dir_2d = lda_dir_2d / norm(lda_dir_2d);
    
    % Convert back to 3D
    second_dir = lda_dir_2d(1) * v1 + lda_dir_2d(2) * v2;
    second_dir = second_dir / norm(second_dir);
    
    % Calculate d-prime
    upper_proj = upper_coords * second_dir;
    lower_proj = lower_coords * second_dir;
    second_dprime = abs(mean(upper_proj) - mean(lower_proj)) / sqrt(0.5 * (var(upper_proj) + var(lower_proj)));
end

%% plot heatmap along three axes
col1=magma(300);
col1=flipud(col1);
f=figure(1);set(f,'Position', [73 422.6000 992 148.8000]);
% til=tiledlayout(length(par_seq),8, 'TileSpacing','none', 'Padding','tight');
% set(til, 'OuterPosition', [0 0 0.8 1]);
hm=heatmap(numbers, 'Colormap',col1,'ColorLimits',[0,1],...%[round(min(min(numbers)),2) round(max(max(numbers)),2)],...
    'GridVisible','off', 'CellLabelColor','none');
hm.Title= 'Separation of medial and lateral input distributions along three main brain axes';
hm.XDisplayLabels = all_labels(:,1);
hm.YDisplayLabels={'ML', 'AP', 'DV'};
s = struct(hm);
s.XAxis.TickLabelRotation = 90; 
set(hm, 'OuterPosition', [0 0 1 0.9]);
tl_pos=get(s.Axes.Title, 'Position');
set(s.Axes.Title, 'FontSize', 11, 'Position', [tl_pos(1) tl_pos(2)-0.2 tl_pos(3)]);
cl=s.Colorbar;
set(cl, 'YDir', 'reverse')
cl.Ticks=[0 0.5 1];
cl.YTickLabel={'1', '0.5', '0'};
cl.Label.String=["separation"; "fraction"];
cl.Label.FontSize=9;


%hm.XDisplayLabels.TickLabelRotation=90;
% hm.YDisplayLabels = nan(size(hm.YDisplayData));
% hm.XDisplayLabels = nan(size(hm.XDisplayData));
% hm.YDisplayLabels = nan(size(hm.YDisplayData));

%% plot heatmap of max value
% Combine data from both celltypes
% numbers_combined = [numbers_all{1}; numbers_all{2}]; % 2 x num_areas matrix

numbers_combined = max_separation_mean; % 2 x num_areas matrix
numbers_combined(cellfun(@isempty, numbers_combined)) = {NaN};
numbers_combined=cell2mat(numbers_combined);

% Find columns where BOTH celltypes have NaN
both_nan = all(isnan(numbers_combined), 1);

% Filter out those columns
numbers_filtered = numbers_combined(:, ~both_nan);
areas_filtered = union_labels(~both_nan, 1);
regions_filtered = union_labels(~both_nan, 3);  % Brain regions (as numbers in cell array)

% Convert cell array of numbers to numeric array
regions_numeric = cell2mat(regions_filtered);

% ===== TWO-LEVEL SORTING =====
% Get unique regions
unique_regions = unique(regions_numeric);

% Get intensity values (first row = celltype 1)
intensity_values = numbers_filtered(1, :);

% Sort within each region by intensity (ascending)
sorted_indices = [];
for r = 1:length(unique_regions)
    in_region = find(regions_numeric == unique_regions(r));
    region_intensities = intensity_values(in_region);
    [~, intensity_order] = sort(region_intensities, 'descend', 'MissingPlacement', 'last');
    sorted_indices = [sorted_indices, in_region(intensity_order)'];
end

% Apply sorting
numbers_sorted = numbers_filtered(:, sorted_indices);
areas_sorted = areas_filtered(sorted_indices);
regions_sorted = regions_numeric(sorted_indices);

% ===== FIND REGION BOUNDARIES =====
region_boundaries = find(diff(regions_sorted) ~= 0);

% Colormap
col1 = magma(300);
% col1 = flipud(col1);

% Create figure
f = figure(1);
set(f, 'Position', [73 441 992 180]);

% Create heatmap with sorted data
hm = heatmap(numbers_sorted, 'Colormap', col1, ...
    'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 10);
hm.Title = 'separation of medial and lateral input distributions';
hm.MissingDataColor = [0.5 0.5 0.5];
hm.XDisplayLabels = areas_sorted;
hm.YDisplayLabels = celltypes;

% Formatting
s = struct(hm);
s.Axes.Title.FontWeight = 'normal';
s.XAxis.TickLabelRotation = 90;
set(hm, 'OuterPosition', [0 0 1 0.8]);
set(s.Axes.Title, 'FontSize', 12);

% Colorbar formatting
cl = s.Colorbar;
% set(cl, 'YDir', 'reverse');
min_val=min(numbers_combined,[],'all');
max_val=max(numbers_combined,[],'all');
med_val=mean([min_val, max_val]);
cl.Ticks = [min_val med_val max_val];
% cl.YTickLabel = {sprintf('%.2f', 1-min_val), sprintf('%.2f', med_val), sprintf('%.2f', 1-max_val)};
cl.YTickLabel = {sprintf('%.2f', min_val), sprintf('%.2f', med_val), sprintf('%.2f', max_val)};
cl.Label.String = ["separation"; "fraction"];
cl.Label.FontSize = 9;

% ===== ADD REGION SEPARATOR LINES USING ANNOTATION =====
% Get heatmap position in figure coordinates
drawnow;  % Force rendering to get correct positions
hm_pos = s.Axes.Position;

num_cols = size(numbers_sorted, 2);
num_rows = size(numbers_sorted, 1);

% Line settings
line_color = [0 0 0];
line_width = 2.5;

for b = 1:length(region_boundaries)
    % Calculate x position in normalized figure coordinates
    x_frac = (region_boundaries(b)) / num_cols;  % fraction along heatmap
    x_fig = hm_pos(1) + x_frac * hm_pos(3);  % convert to figure coordinates
    
    % Draw line from bottom to top of heatmap
    y_bottom = hm_pos(2);
    y_top = hm_pos(2) + hm_pos(4);
    
    annotation(f, 'line', [x_fig x_fig], [y_bottom y_top], ...
        'Color', line_color, 'LineWidth', line_width);
end


% ===== CORRELATION ANALYSIS BETWEEN CELL TYPES =====
celltype1_data = 1-numbers_sorted(1, :);
celltype2_data = 1-numbers_sorted(2, :);

valid_pairs = ~isnan(celltype1_data) & ~isnan(celltype2_data);
ct1_valid = celltype1_data(valid_pairs)';
ct2_valid = celltype2_data(valid_pairs)';
areas_valid = areas_sorted(valid_pairs);

n_valid = sum(valid_pairs);
fprintf('Number of areas with data in both celltypes: %d\n\n', n_valid);

% ===== PEARSON CORRELATION =====
[r_pearson, p_pearson] = corr(ct1_valid, ct2_valid, 'Type', 'Pearson');
fprintf('Pearson correlation:\n');
fprintf('  r = %.3f\n', r_pearson);
fprintf('  p = %.4f\n\n', p_pearson);

% ===== BOOTSTRAP CI FOR PEARSON R =====
nboot = 10000;
boot_r = zeros(nboot, 1);
x_fit = linspace(min(ct1_valid), max(ct1_valid), 100);
y_boot = zeros(nboot, length(x_fit));

rng(42);
for b = 1:nboot
    idx = randi(n_valid, n_valid, 1);  % resample areas
    bx = ct1_valid(idx);
    by = ct2_valid(idx);
    if numel(unique(bx)) > 1 && numel(unique(by)) > 1
        boot_r(b) = corr(bx, by, 'type', 'Pearson');
        p_b = polyfit(bx, by, 1);
        y_boot(b, :) = polyval(p_b, x_fit);
    else
        boot_r(b) = NaN;
        y_boot(b, :) = NaN;
    end
end

ci_lower = prctile(boot_r, 2.5);
ci_upper = prctile(boot_r, 97.5);

y_lower = prctile(y_boot, 2.5, 1);
y_upper = prctile(y_boot, 97.5, 1);
y_fit_mean = median(y_boot, 1, 'omitnan');

fprintf('Bootstrap 95%% CI: [%.3f, %.3f]\n', ci_lower, ci_upper);
fprintf('Bootstrap mean r: %.3f\n', mean(boot_r, 'omitnan'));
fprintf('Bootstrap std: %.3f\n\n', std(boot_r, 'omitnan'));

% ===== SPEARMAN CORRELATION =====
[r_spearman, p_spearman] = corr(ct1_valid, ct2_valid, 'Type', 'Spearman');
fprintf('Spearman correlation:\n');
fprintf('  rho = %.3f\n', r_spearman);
fprintf('  p = %.4f\n\n', p_spearman);

% ===== R-SQUARED =====
r_squared = r_pearson^2;
fprintf('R-squared (variance explained): %.3f (%.1f%%)\n\n', r_squared, r_squared * 100);

% ===== LINEAR REGRESSION =====
coeffs = polyfit(ct1_valid, ct2_valid, 1);
slope = coeffs(1);
intercept = coeffs(2);
fprintf('Linear regression:\n');
fprintf('  %s = %.3f * %s + %.3f\n\n', celltypes{2}, slope, celltypes{1}, intercept);

% ===== PLOT =====
figure; set(gcf, 'position', [329.8000 326.6000 482.4000 344.8000]);
hold on;

% CI band
fill([x_fit, fliplr(x_fit)], [y_lower, fliplr(y_upper)], ...
    [0.3 0.3 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% Regression line
plot(x_fit, y_fit_mean, 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

% Scatter
scatter(ct1_valid, ct2_valid, 30, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 1);

% Labels
label_points = false;
if label_points
    for i = 1:length(ct1_valid)
        text(ct1_valid(i) + 0.01, ct2_valid(i) + 0.01, areas_valid{i}, ...
            'FontSize', 7, 'Color', [0.3 0.3 0.3]);
    end
end

set(gca, 'FontSize', 12)
xlabel(sprintf('%s (separation fraction)', celltypes{1}));
ylabel(sprintf('%s (separation fraction)', celltypes{2}));

tit = title(sprintf('r=%.2f [%.2f, %.2f], p=%.3g', r_pearson, ci_lower, ci_upper, p_pearson));
tit.FontWeight = 'normal';
tit.FontSize = 12;

axis square;
hold off;

% ===== BOOTSTRAP DISTRIBUTION =====
figure;
histogram(boot_r, 50, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.6]);
hold on;
xline(r_pearson, 'r-', 'LineWidth', 2, 'Label', sprintf('r=%.2f', r_pearson));
xline(ci_lower, 'b--', 'LineWidth', 1.5);
xline(ci_upper, 'b--', 'LineWidth', 1.5);
xlabel('Pearson r');
ylabel('Probability');
title(sprintf('Bootstrap Distribution: r=%.2f [%.2f, %.2f]', r_pearson, ci_lower, ci_upper));
set(gca, 'FontSize', 12);

% ===== SUMMARY =====
correlation_results = table(...
    {'Pearson'; 'Spearman'}, ...
    [r_pearson; r_spearman], ...
    [p_pearson; p_spearman], ...
    'VariableNames', {'Method', 'Correlation', 'P_value'});
disp('Summary:');
disp(correlation_results);


%% plot 2D scatter of single areas with ksdensity distributions
% prepare data

dt_all = {};
myCoord_all = {};
for ct = 1:length(celltypes)
    dt={};
    for t = 1:length(myType)
        typeName = myType{t};
        dt{t} = load(fullfile(spath,[celltypes{ct}, '_',typeName,'.mat'])); 
    end
    dt_all{ct} = dt;

    labels = dt{1}.TT(:,1:3);
    myCoord={};
    min_exp_nr=min(length(dt{1,1}.expid), length(dt{1,2}.expid));

    if length(dt{1,1}.expid) == min_exp_nr
        for i=1:min_exp_nr
            id=dt{1,1}.expid{i}(1:9);
            index_2=find(contains(dt{1,2}.expid, id));
            if index_2
                myCoord{1,1}(:,i)=dt{1,1}.CoordMerge(:, i);
                myCoord{1,2}(:,i)=dt{1,2}.CoordMerge(:, index_2);
            end
        end
    else
        for i=1:min_exp_nr
            id=dt{1,2}.expid{i}(1:9);
            index_1=find(contains(dt{1,1}.expid, id));
            if index_1
                myCoord{1,1}(:,i)=dt{1,1}.CoordMerge(:, index_1);
                myCoord{1,2}(:,i)=dt{1,2}.CoordMerge(:, i);
            end
        end
    end
    myCoord_all{ct} = myCoord;
end
t_per_ct = [3, 1];

%% plot them all within one figure
area='ZI';
AX1=3;
AX2=2;

for ct = 1:length(celltypes)
    t = t_per_ct(ct);
    labels = dt_all{ct}{1}.TT(:,1:3);
    myCoord = myCoord_all{ct};

    Coordinates_tot=[];
    coord_tot={};

    for tt = 1:length(myType)
        typeName = myType{tt};
        idSC = find(ismember(labels(:,1),area));
        myCoord_temp=myCoord{tt}(idSC,:);

        [sa, sl] = size(myCoord_temp);coordXYZ = [];
        Coord_type={};
        for i = 1:sl
            coordXYZ = cell2mat(myCoord_temp(:,i));
            if ~isempty(coordXYZ)
                coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
                coordXYZ=coordXYZ+randn(size(coordXYZ))/300;
                coordXYZ(:,1)=coordXYZ(:,1)+(-0.05 + (0.05+0.05).*rand(size(coordXYZ,1),1));
                Coordinates_tot=[Coordinates_tot; coordXYZ];
                Coord_type{i}=coordXYZ;
            else
                coordXYZ=[NaN NaN NaN];
                Coord_type{i}=coordXYZ;
            end
        end
        coord_tot{1,tt}=Coord_type;
    end

    f=figure();
    bigT=tiledlayout('flow','TileSpacing','compact', 'Padding','compact');
    bigT.Title.String = [area ' - ' celltypes{ct}];
    bigT.Title.FontSize = 14;
    bigT.Title.FontWeight = 'bold';
    counter=1;
    for fi=1:length(coord_tot{1})
        if size(coord_tot{1,1}{fi},1)>=20 && size(coord_tot{1,2}{fi},1)>=20
            nexttile(bigT); set(gca, 'box', 'off', 'visible','off')
            fact=2;
            til=tiledlayout(bigT,3*fact,3*fact, 'TileSpacing','none', 'Padding','tight');
            til.Layout.Tile = counter;
            til.Title.String=area;
            til.Title.FontSize=14;
            til.Title.FontWeight='bold';

            nexttile(til,[1 2*fact+1]);
            [y_med,~] = ksdensity(coord_tot{1,1}{fi}(:,AX1),'Support', [min(coord_tot{1,1}{fi}(:,AX1))-0.1 ...
                max(coord_tot{1,1}{fi}(:,AX1))+0.1], 'NumPoints',300, 'BandWidth', 0.3);
            x_med=linspace(min(coord_tot{1,1}{fi}(:,AX1)), max(coord_tot{1,1}{fi}(:,AX1)),300);
            plot(x_med,y_med,'LineWidth',2, 'Color',col(t,:)); hold on;
            [y_lat,~] = ksdensity(coord_tot{1,2}{fi}(:,AX1), 'Support', [min(coord_tot{1,2}{fi}(:,AX1))-0.1 ...
                max(coord_tot{1,2}{fi}(:,AX1))+0.1], 'NumPoints',300,'BandWidth', 0.3);
            x_lat=linspace(min(coord_tot{1,2}{fi}(:,AX1)), max(coord_tot{1,2}{fi}(:,AX1)),300);
            plot(x_lat,y_lat,'LineWidth',2, 'Color',col(t+1,:));

            newx=linspace(min([min(coord_tot{1,1}{fi}(:,AX1)) min(coord_tot{1,2}{fi}(:,AX1))]), ...
                max([max(coord_tot{1,1}{fi}(:,AX1)) max(coord_tot{1,2}{fi}(:,AX1))]),300);
            new_y_med=interp1(x_med,y_med,newx);
            new_y_lat=interp1(x_lat,y_lat,newx);
            pmin=min(new_y_med,new_y_lat);
            pmax=max(new_y_med,new_y_lat);

            overlap_horizontal=sum(pmin)./sum(pmax);
            [~, pval_horizontal,stat_horizontal] = kstest2(coord_tot{1,1}{fi}(:,AX1), coord_tot{1,2}{fi}(:,AX1));

            xlim(gca,[(min(Coordinates_tot(:,AX1))-0.2) (max(Coordinates_tot(:,AX1))+0.2)]);
            set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');

            nexttile(til,4*fact-1,[2*fact+1 2*fact+1]);
            s_med=scatter(coord_tot{1,1}{fi}(:,AX1), coord_tot{1,1}{fi}(:,AX2), 10, 'filled', CData=col(t,:), MarkerFaceAlpha=0.6); hold on;
            s_lat=scatter(coord_tot{1,2}{fi}(:,AX1), coord_tot{1,2}{fi}(:,AX2), 10, 'filled', CData=col(t+1,:), MarkerFaceAlpha=0.6);hold on;
            uistack(s_lat, 'bottom');
            set(gca,'box', 'on', 'Color', 'none');
            xlim(gca,[(min(Coordinates_tot(:,AX1))-0.2) (max(Coordinates_tot(:,AX1))+0.2)]);
            ylim(gca,[(min(Coordinates_tot(:,AX2))-0.2) (max(Coordinates_tot(:,AX2))+0.2)]);
            if AX2==2
                set(gca, 'YDir', 'reverse')
            end
            ax=gca; ax.FontSize=10;

            if AX1==1
                xaxis='A-P';
                ax.XTickLabel = split(num2str(xticks*(-1)))';
            elseif AX1==2
                xaxis='D-V';
                ax.XTickLabel = split(num2str(xticks*(-1)))';
            else
                xaxis='M-L';
                ax.XTickLabel = split(num2str(xticks))';
            end

            if AX2==1
                yaxis='A-P';
                ax.YTickLabel = split(num2str(yticks*(-1)))';
            elseif AX2==2
                yaxis='D-V';
                ax.YTickLabel = split(num2str(yticks*(-1)))';
            else
                yaxis='M-L';
                ax.YTickLabel = split(num2str(yticks))';
            end

            xlabel(xaxis, 'FontSize', 14); ylabel(yaxis,'FontSize', 14);
            set(gca, 'FontSize', 12);

            nexttile(til,6*fact, [2*fact+1 1])
            [y_med_vertical,~] = ksdensity(coord_tot{1,1}{fi}(:,AX2),'Support', [min(coord_tot{1,1}{fi}(:,AX2))-0.1 ...
                max(coord_tot{1,1}{fi}(:,AX2))+0.1], 'NumPoints',300, 'BandWidth', 0.3);
            x_med_vertical=linspace(min(coord_tot{1,1}{fi}(:,AX2)), max(coord_tot{1,1}{fi}(:,AX2)),300);
            plot(y_med_vertical,x_med_vertical,'LineWidth',2, 'Color',col(t,:)); hold on;
            [y_lat_vertical,~] = ksdensity(coord_tot{1,2}{fi}(:,AX2),'Support', [min(coord_tot{1,2}{fi}(:,AX2))-0.1 ...
                max(coord_tot{1,2}{fi}(:,AX2))+0.1], 'NumPoints',300,'BandWidth', 0.3);
            x_lat_vertical=linspace(min(coord_tot{1,2}{fi}(:,AX2)), max(coord_tot{1,2}{fi}(:,AX2)),300);
            plot(y_lat_vertical,x_lat_vertical,'LineWidth',2, 'Color',col(t+1,:));

            newx_vertical=linspace(min([min(coord_tot{1,1}{fi}(:,AX2)) min(coord_tot{1,2}{fi}(:,AX2))]), ...
                max([max(coord_tot{1,1}{fi}(:,AX2)) max(coord_tot{1,2}{fi}(:,AX2))]),300);
            new_y_med_vertical=interp1(x_med_vertical,y_med_vertical,newx_vertical);
            new_y_lat_vertical=interp1(x_lat_vertical,y_lat_vertical,newx_vertical);
            pmin_vertical=min(new_y_med_vertical,new_y_lat_vertical);
            pmax_vertical=max(new_y_med_vertical,new_y_lat_vertical);

            overlap_vertical=sum(pmin_vertical)./sum(pmax_vertical);
            [~, pval_vertical,stat_vertical] = kstest2(coord_tot{1,1}{fi}(:,AX2), coord_tot{1,2}{fi}(:,AX2));

            ylim(gca,[(min(Coordinates_tot(:,AX2))-0.2) (max(Coordinates_tot(:,AX2))+0.2)]);
            set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
            if AX2==2
                set(gca, 'YDir', 'reverse')
            end

            x_range=(max(Coordinates_tot(:,AX1))+0.2) -(min(Coordinates_tot(:,AX1))-0.2);
            y_range=(max(Coordinates_tot(:,AX2))+0.2) -(min(Coordinates_tot(:,AX2))-0.2);
            xy_ratio=y_range./x_range;

            counter=counter+1;
        end
    end
    fig_width=100*counter;
    set(f, "Position", [29 150 fig_width fig_width*xy_ratio], 'renderer','Painters');
end

%% plot them all within one graph

area='ZI';
AX1=3;
AX2=2;

for ct = 1:length(celltypes)
    dt={};
    for tt = 1:length(myType)
        typeName = myType{tt};
        dt{tt} = load(fullfile(spath,[celltypes{ct}, '_',typeName,'.mat'])); 
    end
    t = t_per_ct(ct);
    labels = dt{1}.TT(:,1:3);
    myCoord={};
    myCoord{1,1}(:,:)=dt{1,1}.CoordMerge(:,:);
    myCoord{1,2}(:,:)=dt{1,2}.CoordMerge(:,:);

    Coordinates_tot=[];
    Coordinates_tot_type=cell(1, length(myType));
    coord_tot={};

    for tt = 1:length(myType)
        coll = col(tt+2,:);
        typeName = myType{tt};
        idSC = find(ismember(labels(:,1),area));
        myCoord_temp=myCoord{tt}(idSC,:);

        [sa, sl] = size(myCoord_temp);coordXYZ = [];
        Coord_type={};
        for i = 1:sl
            coordXYZ = cell2mat(myCoord_temp(:,i));
            if ~isempty(coordXYZ)
                coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
                coordXYZ=coordXYZ+randn(size(coordXYZ))/300;
                coordXYZ(:,1)=coordXYZ(:,1)+(-0.05 + (0.05+0.05).*rand(size(coordXYZ,1),1));
                Coordinates_tot=[Coordinates_tot; coordXYZ];
                Coordinates_tot_type{tt}=[Coordinates_tot_type{tt};coordXYZ];
                Coord_type{i}=coordXYZ;
            else
                coordXYZ=[NaN NaN NaN];
                Coord_type{i}=coordXYZ;
            end
        end
        coord_tot{1,tt}=Coord_type;
    end

    f=figure();
    hold on;

    fact=2;
    til=tiledlayout(3*fact,3*fact, 'TileSpacing','none', 'Padding','tight');
    til.Title.String=[area ' - ' celltypes{ct}];
    til.Title.FontSize=14;
    til.Title.FontWeight='bold';

    nexttile(til,[1 2*fact+1]);
    [y_med,~] = ksdensity(Coordinates_tot_type{1}(:,AX1),'Support', [min(Coordinates_tot_type{1}(:,AX1))-0.1 ...
        max(Coordinates_tot_type{1}(:,AX1))+0.1], 'NumPoints',300, 'BandWidth', 0.25);
    x_med=linspace(min(Coordinates_tot_type{1}(:,AX1)), max(Coordinates_tot_type{1}(:,AX1)),300);
    plot(x_med,y_med,'LineWidth',2, 'Color',col(t,:)); hold on;
    [y_lat,~] = ksdensity(Coordinates_tot_type{2}(:,AX1), 'Support', [min(Coordinates_tot_type{2}(:,AX1))-0.1 ...
        max(Coordinates_tot_type{2}(:,AX1))+0.1], 'NumPoints',300,'BandWidth', 0.25);
    x_lat=linspace(min(Coordinates_tot_type{2}(:,AX1)), max(Coordinates_tot_type{2}(:,AX1)),300);
    plot(x_lat,y_lat,'LineWidth',2, 'Color',col(t+1,:));

    newx=linspace(min([min(Coordinates_tot_type{1}(:,AX1)) min(Coordinates_tot_type{2}(:,AX1))]), ...
        max([max(Coordinates_tot_type{1}(:,AX1)) max(Coordinates_tot_type{2}(:,AX1))]),300);
    new_y_med=interp1(x_med,y_med,newx);
    new_y_lat=interp1(x_lat,y_lat,newx);
    pmin=min(new_y_med,new_y_lat);
    pmax=max(new_y_med,new_y_lat);

    overlap_horizontal=sum(pmin)./sum(pmax);
    [~, pval_horizontal,stat_horizontal] = kstest2(Coordinates_tot_type{1}(:,AX1), Coordinates_tot_type{2}(:,AX1));

    xlim(gca,[(min(Coordinates_tot(:,AX1))-0.2) (max(Coordinates_tot(:,AX1))+0.2)]);
    set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');

    nexttile(til,4*fact-1,[2*fact+1 2*fact+1]);
    s_med=scatter(Coordinates_tot_type{1}(:,AX1), Coordinates_tot_type{1}(:,AX2), 10, 'filled', CData=col(t,:), MarkerFaceAlpha=0.6); hold on;
    s_lat=scatter(Coordinates_tot_type{2}(:,AX1), Coordinates_tot_type{2}(:,AX2), 10, 'filled', CData=col(t+1,:), MarkerFaceAlpha=0.6);hold on;
    uistack(s_lat, 'bottom');

    if AX2==2
        set(gca, 'YDir', 'reverse')
    end
    ax=gca; ax.FontSize=10;
    xlim(ax,[(min(Coordinates_tot(:,AX1))-0.2) (max(Coordinates_tot(:,AX1))+0.2)]);
    ylim(ax,[(min(Coordinates_tot(:,AX2))-0.2) (max(Coordinates_tot(:,AX2))+0.2)]);
    set(ax,'box', 'on', 'Color', 'none');

    if AX1==1
        xaxis='A-P';
        ax.XTickLabel = split(num2str(xticks*(-1)))';
    elseif AX1==2
        xaxis='D-V';
        ax.XTickLabel = split(num2str(xticks*(-1)))';
    else
        xaxis='M-L';
        ax.XTickLabel = split(num2str(xticks))';
    end

    if AX2==1
        yaxis='A-P';
        ax.YTickLabel = split(num2str(yticks*(-1)))';
    elseif AX2==2
        yaxis='D-V';
        ax.YTickLabel = split(num2str(yticks*(-1)))';
    else
        yaxis='M-L';
        ax.YTickLabel = split(num2str(yticks))';
    end

    xlabel(xaxis, 'FontSize', 14); ylabel(yaxis,'FontSize', 14);
    set(gca, 'FontSize', 12);

    nexttile(til,6*fact, [2*fact+1 1])
    [y_med_vertical,~] = ksdensity(Coordinates_tot_type{1}(:,AX2),'Support', [min(Coordinates_tot_type{1}(:,AX2))-0.1 ...
        max(Coordinates_tot_type{1}(:,AX2))+0.1], 'NumPoints',300, 'BandWidth', 0.25);
    x_med_vertical=linspace(min(Coordinates_tot_type{1}(:,AX2)), max(Coordinates_tot_type{1}(:,AX2)),300);
    plot(y_med_vertical,x_med_vertical,'LineWidth',2, 'Color',col(t,:)); hold on;
    [y_lat_vertical,~] = ksdensity(Coordinates_tot_type{2}(:,AX2),'Support', [min(Coordinates_tot_type{2}(:,AX2))-0.1 ...
        max(Coordinates_tot_type{2}(:,AX2))+0.1], 'NumPoints',300,'BandWidth', 0.25);
    x_lat_vertical=linspace(min(Coordinates_tot_type{2}(:,AX2)), max(Coordinates_tot_type{2}(:,AX2)),300);
    plot(y_lat_vertical,x_lat_vertical,'LineWidth',2, 'Color',col(t+1,:));

    newx_vertical=linspace(min([min(Coordinates_tot_type{1}(:,AX2)) min(Coordinates_tot_type{2}(:,AX2))]), ...
        max([max(Coordinates_tot_type{1}(:,AX2)) max(Coordinates_tot_type{2}(:,AX2))]),300);
    new_y_med_vertical=interp1(x_med_vertical,y_med_vertical,newx_vertical);
    new_y_lat_vertical=interp1(x_lat_vertical,y_lat_vertical,newx_vertical);
    pmin_vertical=min(new_y_med_vertical,new_y_lat_vertical);
    pmax_vertical=max(new_y_med_vertical,new_y_lat_vertical);

    overlap_vertical=sum(pmin_vertical)./sum(pmax_vertical);
    [~, pval_vertical,stat_vertical] = kstest2(Coordinates_tot_type{1}(:,AX2), Coordinates_tot_type{2}(:,AX2));

    ylim(gca,[(min(Coordinates_tot(:,AX2))-0.2) (max(Coordinates_tot(:,AX2))+0.2)]);
    set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
    if AX2==2
        set(gca, 'YDir', 'reverse')
    end
    x_range=(max(Coordinates_tot(:,AX1))+0.3) -(min(Coordinates_tot(:,AX1))-0.3);
    y_range=(max(Coordinates_tot(:,AX2))+0.3) -(min(Coordinates_tot(:,AX2))-0.3);
    xy_ratio=y_range./x_range;
    fig_width=280;

    set(f, "Position", [29 150 fig_width fig_width*xy_ratio], 'renderer','Painters');
end

%% plot them all within one graph - MED LAT AS ONE

dt={};
for t = 1:length(myType)
    typeName = myType{t};
    dt{t} = load(fullfile(spath,[celltypes{1}, '_',typeName,'.mat'])); 
end

labels = dt{1}.TT(:,1:3);
myCoord={};

myCoord{1,1}(:,:)=dt{1,1}.CoordMerge(:,:);
myCoord{1,2}(:,:)=dt{1,2}.CoordMerge(:,:);


Coordinates_tot=[];
Coordinates_tot_type=cell(1, length(myType));
coord_tot={};

area='PAG';

AX1=3; %distrubition along x-axis; ML=3; AP=1; DV=2;
AX2=2; %distribution along y-axis

for t = 1:length(myType)
    typeName = myType{t};
    idSC = find(ismember(labels(:,1),area));
    myCoord_temp=myCoord{t}(idSC,:);
%     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
%         if threCur(i,1)==0
%             myCoord(i,:)={[]};
%         end
%     end

    %plot transformed coordinates
    %figure(f);
    [sa, sl] = size(myCoord_temp);coordXYZ = [];
    %     theta=deg2rad(-28.022);
    %     phi=deg2rad(-8.29);
    Coord_type={};
    for i = 1:sl
        coordXYZ = cell2mat(myCoord_temp(:,i));
        if ~isempty(coordXYZ)
            coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:); %keep only ipsilateral
            coordXYZ=coordXYZ+randn(size(coordXYZ))/300;%add some small random noise to visualize scatter points that would be on the same pixel
            coordXYZ(:,1)=coordXYZ(:,1)+(-0.05 + (0.05+0.05).*rand(size(coordXYZ,1),1)); %jitter points along AP axis to cover slice spacing effect
            %         X=coordXYZ(:,1);
            %         Y=coordXYZ(:,3);
            %         Z=coordXYZ(:,2);
            %         Xnew=X*cos(phi)-Y*sin(theta)*sin(phi)-Z*cos(theta)*sin(phi);
            %         Ynew=Y*cos(theta)-Z*sin(theta);
            %         Znew=X*sin(phi)+Y*sin(theta)*cos(phi)+Z*cos(theta)*cos(phi);
            %         coordXYZ=[Xnew Znew Ynew];
            Coordinates_tot=[Coordinates_tot; coordXYZ];
            Coordinates_tot_type{t}=[Coordinates_tot_type{t};coordXYZ];
            Coord_type{i}=coordXYZ;
        else
        coordXYZ=[NaN NaN NaN];
        Coord_type{i}=coordXYZ;
        end
    end
    coord_tot{1,t}=Coord_type;
end

t=5; %5 for ntsr, 6 for grp

f=figure();
hold on;
fact=2;
til=tiledlayout(3*fact,3*fact, 'TileSpacing','none', 'Padding','tight');
set(til, 'OuterPosition', [0 0 0.98 0.98]);
til.Title.String=area;
til.Title.FontSize=14;
til.Title.FontWeight='bold';
%plot horizontal density curve
nexttile(til,[1 2*fact+1]);

[y_med,~] = ksdensity(Coordinates_tot(:,AX1),'Support', [min(Coordinates_tot(:,AX1))-0.1 ...
    max(Coordinates_tot(:,AX1))+0.1], 'NumPoints',300, 'BandWidth', 0.2); %ML
x_med=linspace(min(Coordinates_tot(:,AX1)), max(Coordinates_tot(:,AX1)),300);
plot(x_med,y_med,'LineWidth',2, 'Color',col(t,:)); hold on;


xlim(gca,[(min(Coordinates_tot(:,AX1))-0.2) (max(Coordinates_tot(:,AX1))+0.2)]);
set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');

nexttile(til,4*fact-1,[2*fact+1 2*fact+1]);
hold on;

ax=gca; ax.FontSize=10;
% axis(ax, 'equal');
xlim(ax,[(min(Coordinates_tot(:,AX1))-0.2) (max(Coordinates_tot(:,AX1))+0.2)]);
ylim(ax,[(min(Coordinates_tot(:,AX2))-0.2) (max(Coordinates_tot(:,AX2))+0.2)]);
set(ax,'box', 'on', 'Color', 'none');

%----------------DENSITY MAP INSTEAD OF SCATTER----
% compute 2D density
x = Coordinates_tot(:,AX1);
y = Coordinates_tot(:,AX2);

% define grid (300x300 is fine)
numPoints = 300;
xlin = linspace(min(x)-0.2, max(x)+0.2, numPoints);
ylin = linspace(min(y)-0.2, max(y)+0.2, numPoints);
[Xgrid, Ygrid] = meshgrid(xlin, ylin);

% evaluate 2D kernel density
[ks,~] = ksdensity([x y],[Xgrid(:) Ygrid(:)], ...
    'Bandwidth',[0.2 0.08], ...      % use your previous bandwidths
    'Support','unbounded');          % or specify exact range

F = reshape(ks, numPoints, numPoints);
F(F == 0) = NaN;  

% % plot density map instead of scatter
h=imagesc(xlin, ylin, F); 
set(gca, 'YDir', 'normal'); % so Y increases upward
nColors = 256;           % number of colormap levels
targetColor = col(t,:)-0.2;  % your desired color, e.g., [0 0.4470 0.7410]

% interpolate from white to the color
cmap = [linspace(1,targetColor(1),nColors)', ...
        linspace(1,targetColor(2),nColors)', ...
        linspace(1,targetColor(3),nColors)'];

colormap(cmap);



if AX2==2
    set(gca, 'YDir', 'reverse')
end


%-------------------------------


% %plot scatter points
% s_med=scatter(Coordinates_tot(:,AX1), Coordinates_tot(:,AX2), 3, 'filled', CData=col(t,:), MarkerFaceAlpha=0.3); hold on;

% % optional: overlay contour lines
% contour(Xgrid, Ygrid, F, 8, 'k', 'LineWidth', 0.5);
% if AX2==2
%     set(gca, 'YDir', 'reverse')
% end

if AX1==1
    xaxis='A-P';
    ax.XTickLabel = split(num2str(xticks*(-1)))';
elseif AX1==2
    xaxis='D-V';
    ax.XTickLabel = split(num2str(xticks*(-1)))';
else
    xaxis='M-L';
    ax.XTickLabel = split(num2str(xticks))';
end

if AX2==1
    yaxis='A-P';
    ax.YTickLabel = split(num2str(yticks*(-1)))';
elseif AX2==2
    yaxis='D-V';
    ax.YTickLabel = split(num2str(yticks*(-1)))';
else
    yaxis='M-L';
    ax.YTickLabel = split(num2str(yticks))';
end

xlabel(xaxis, 'FontSize', 14); ylabel(yaxis,'FontSize', 14);
% axis(gca,'equal');
set(gca, 'FontSize', 12);

%plot vertical density curve
nexttile(til,6*fact, [2*fact+1 1])
[y_med_vertical,~] = ksdensity(Coordinates_tot(:,AX2),'Support', [min(Coordinates_tot(:,AX2))-0.1 ...
    max(Coordinates_tot(:,AX2))+0.1],'NumPoints',300, 'BandWidth', 0.08);
x_med_vertical=linspace(min(Coordinates_tot(:,AX2)), max(Coordinates_tot(:,AX2)),300);
plot(y_med_vertical,x_med_vertical,'LineWidth',2, 'Color',col(t,:)); hold on;


ylim(gca,[(min(Coordinates_tot(:,AX2))-0.2) (max(Coordinates_tot(:,AX2))+0.2)]);
set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
if AX2==2
    set(gca, 'YDir', 'reverse')
end
x_range=(max(Coordinates_tot(:,AX1))+0.3) -(min(Coordinates_tot(:,AX1))-0.3);
y_range=(max(Coordinates_tot(:,AX2))+0.3) -(min(Coordinates_tot(:,AX2))-0.3);
xy_ratio=y_range./x_range;
fig_width=214;

set(f, "Position", [29 150 fig_width fig_width*xy_ratio], 'renderer','Painters');

