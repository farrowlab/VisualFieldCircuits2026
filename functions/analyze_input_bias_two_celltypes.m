function Results = analyze_input_bias_two_celltypes(nu1, nl1, nu2, nl2, areaNames, varargin)
%ANALYZE_INPUT_BIAS_TWO_CELLTYPES Upper vs lower input bias of two cell types and their difference.
%   Results = ANALYZE_INPUT_BIAS_TWO_CELLTYPES(nu1, nl1, nu2, nl2, areaNames, Name, Value, ...)
%
%   nu1, nl1   animals-by-areas cell counts of cell type 1 (CT1) for the
%              upper (U, medial sSC) and lower (L, lateral sSC) injections
%   nu2, nl2   the same for cell type 2 (CT2)
%   areaNames  area acronyms, one per column. An area whose counts are NaN
%              for a cell type is skipped for that cell type.
%
%   Per area and cell type the effect is the log odds ratio of the input
%   fraction, U vs L (> 0: bias towards the upper visual field), tested by
%   a bootstrap over animals (BCa CIs) and a permutation test of the U/L
%   labels. The interaction is the difference of the two log odds ratios
%   (CT1 - CT2), tested by a bootstrap and by permuting the cell-type
%   labels within the U and within the L animals.
%
%   Name-value options [default]
%     TotalU1, TotalL1, TotalU2, TotalL2   total inputs per animal [row sums]
%     Weighted          false: mean of per-animal proportions;
%                       true: pooled proportions [false]
%     NBootstrap        bootstrap resamples [10000]
%     NPermutations     permutations [10000]
%     Seed              random seed [42]
%     Alpha             pseudocount for log2FC [1e-10]
%     QDiscovery        q cutoffs counting discoveries for the FCR-adjusted
%                       CIs, [CT1 CT2 interaction] [0.045 0.03 0.045]
%     QColor            volcano plots: q cutoff for colouring, [CT1 CT2 int] [0.05 0.04 0.05]
%     PLabel            volcano plots: p cutoff for labelling, [CT1 CT2 int] [0.052 0.05 0.05]
%     ForestEffectMin   forest plots: |effect| above which a point whose CI
%                       excludes 0 gets its sign colour [0.4]
%     Plot              draw figures [true]
%     ColorU1, ColorL1, ColorU2, ColorL2, ColorIntPos, ColorIntNeg, ColorNonSig
%
%   Results has one row per area. CT q-values and FCR-adjusted CIs are
%   based on the permutation p-values.

% -------------------- Options --------------------
p = inputParser;
addRequired(p, 'nu1', @isnumeric);
addRequired(p, 'nl1', @isnumeric);
addRequired(p, 'nu2', @isnumeric);
addRequired(p, 'nl2', @isnumeric);
addRequired(p, 'areaNames', @(x) iscellstr(x) || isstring(x) || ischar(x));
totalCheck = @(x) isempty(x) || (isnumeric(x) && isvector(x));
addParameter(p, 'TotalU1', [], totalCheck);
addParameter(p, 'TotalL1', [], totalCheck);
addParameter(p, 'TotalU2', [], totalCheck);
addParameter(p, 'TotalL2', [], totalCheck);
addParameter(p, 'Weighted', false, @islogical);
addParameter(p, 'NBootstrap', 10000, @isnumeric);
addParameter(p, 'NPermutations', 10000, @isnumeric);
addParameter(p, 'Seed', 42, @isnumeric);
addParameter(p, 'Alpha', 1e-10, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'QDiscovery', [0.045 0.03 0.045], @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'QColor', [0.05 0.04 0.05], @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'PLabel', [0.052 0.05 0.05], @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'ForestEffectMin', 0.4, @isnumeric);
addParameter(p, 'Plot', true, @islogical);
addParameter(p, 'ColorU1', [0.91 0.68 0.49], @validateColor);
addParameter(p, 'ColorL1', [0.43 0.61 0.77], @validateColor);
addParameter(p, 'ColorU2', [0.91 0.68 0.49], @validateColor);
addParameter(p, 'ColorL2', [0.43 0.61 0.77], @validateColor);
addParameter(p, 'ColorIntPos', [0.91 0.68 0.49], @validateColor);
addParameter(p, 'ColorIntNeg', [0.43 0.61 0.77], @validateColor);
addParameter(p, 'ColorNonSig', [0.78 0.78 0.78], @validateColor);
parse(p, nu1, nl1, nu2, nl2, areaNames, varargin{:});
opt = p.Results;

weighted   = opt.Weighted;
nBootstrap = opt.NBootstrap;
nPerm      = opt.NPermutations;
eps_prop   = 1e-10;            % keeps proportions inside (0, 1)

if weighted
    fprintf('Using WEIGHTED estimation (pooled proportions - biased toward high-trial animals)\n');
else
    fprintf('Using UNWEIGHTED estimation (mean of proportions - unbiased animal-level inference)\n');
end

% -------------------- Sizes, totals, fractions --------------------
areaNames = cellstr(areaNames(:));
nAreas = numel(areaNames);
[nU1, cU1] = size(nu1); [nL1, cL1] = size(nl1);
[nU2, cU2] = size(nu2); [nL2, cL2] = size(nl2);
if ~(cU1 == cL1 && cU1 == cU2 && cU1 == cL2 && cU1 == nAreas)
    error(['All count matrices must have the same number of columns (areas). ' ...
        'areaNames: %d, columns: %d %d %d %d'], nAreas, cU1, cL1, cU2, cL2);
end

N_U1 = group_totals(opt.TotalU1, nu1, nU1, 'TotalU1');
N_L1 = group_totals(opt.TotalL1, nl1, nL1, 'TotalL1');
N_U2 = group_totals(opt.TotalU2, nu2, nU2, 'TotalU2');
N_L2 = group_totals(opt.TotalL2, nl2, nL2, 'TotalL2');
if any(N_U1 == 0) || any(N_L1 == 0) || any(N_U2 == 0) || any(N_L2 == 0)
    error('TotalU/TotalL must be > 0 for all animals.');
end

meanFracU_CT1 = mean(nu1 ./ N_U1, 1)';
meanFracL_CT1 = mean(nl1 ./ N_L1, 1)';
meanFracU_CT2 = mean(nu2 ./ N_U2, 1)';
meanFracL_CT2 = mean(nl2 ./ N_L2, 1)';
log2FC_CT1 = log2((meanFracU_CT1 + opt.Alpha) ./ (meanFracL_CT1 + opt.Alpha));
log2FC_CT2 = log2((meanFracU_CT2 + opt.Alpha) ./ (meanFracL_CT2 + opt.Alpha));

% -------------------- Per-area bootstrap and permutation tests --------------------
% Random numbers are drawn in a fixed order per area: CT1 (bootstrap, then
% permutations), CT2, interaction.
S = [empty_result(nAreas), empty_result(nAreas), empty_result(nAreas)];   % CT1, CT2, interaction
rng(opt.Seed);

for i = 1:nAreas
    s1U = nu1(:, i); s1L = nl1(:, i);
    s2U = nu2(:, i); s2L = nl2(:, i);
    valid1 = ~any(isnan(s1U)) && ~any(isnan(s1L)) && ~(sum(s1U) == 0 && sum(s1L) == 0);
    valid2 = ~any(isnan(s2U)) && ~any(isnan(s2L)) && ~(sum(s2U) == 0 && sum(s2L) == 0);

    if valid1
        try
            S(1) = store_result(S(1), i, bias_boot_perm(s1U, s1L, N_U1, N_L1, nBootstrap, nPerm, weighted, eps_prop));
        catch ME
            warning('CT1 bootstrap/permutation failed for area %d: %s', i, ME.message);
        end
    end
    if valid2
        try
            S(2) = store_result(S(2), i, bias_boot_perm(s2U, s2L, N_U2, N_L2, nBootstrap, nPerm, weighted, eps_prop));
        catch ME
            warning('CT2 bootstrap/permutation failed for area %d: %s', i, ME.message);
        end
    end
    if valid1 && valid2
        try
            S(3) = store_result(S(3), i, interaction_boot_perm(s1U, s1L, s2U, s2L, ...
                N_U1, N_L1, N_U2, N_L2, nBootstrap, nPerm, weighted, eps_prop));
        catch ME
            warning('Interaction bootstrap/permutation failed for area %d: %s', i, ME.message);
        end
    end
end

% -------------------- q-values and adjusted CIs --------------------
% FCR-adjusted CIs: alpha = 0.05 * R / m, with R discoveries at
% q < QDiscovery (Bonferroni if R = 0). Bonferroni: alpha = 0.05 / nAreas.
setNames = {'CT1', 'CT2', 'Interaction'};
R_disc = zeros(1, 3);
for k = 1:3
    S(k).perm_q = bh_qvalues(S(k).perm_p);
    valid = ~isnan(S(k).perm_p);
    m = sum(valid);
    R_disc(k) = sum(S(k).perm_q(valid) < opt.QDiscovery(k));
    if R_disc(k) == 0
        alpha_adj = 0.05 / m;
        fprintf('%s FCR: No discoveries. Using Bonferroni.\n', setNames{k});
    else
        alpha_adj = 0.05 * (R_disc(k) / m);
        fprintf('%s FCR: %d discoveries. Adjusted alpha = %.5f.\n', setNames{k}, R_disc(k), alpha_adj);
    end

    S(k).CIlo_fcr = nan(nAreas, 1);  S(k).CIhi_fcr = nan(nAreas, 1);
    S(k).CIlo_bonf = nan(nAreas, 1); S(k).CIhi_bonf = nan(nAreas, 1);
    for i = 1:nAreas
        if ~isempty(S(k).dist{i}) && ~isnan(S(k).z0(i))
            [S(k).CIlo_fcr(i), S(k).CIhi_fcr(i)] = bca_ci(S(k).dist{i}, S(k).z0(i), S(k).a(i), alpha_adj);
            [S(k).CIlo_bonf(i), S(k).CIhi_bonf(i)] = bca_ci(S(k).dist{i}, S(k).z0(i), S(k).a(i), 0.05 / nAreas);
        end
    end
end

% -------------------- Summary --------------------
excl0 = @(lo, hi) sum((lo > 0) | (hi < 0));
fprintf('\n=======================================================\n');
fprintf('BOOTSTRAP + PERMUTATION SUMMARY\n');
fprintf('=======================================================\n');
for k = 1:3
    fprintf('%s:\n', setNames{k});
    fprintf('  Permutation FDR discoveries: %d\n', R_disc(k));
    fprintf('  Bootstrap BCa (unadj) CI excl. 0: %d\n', excl0(S(k).CIlo_bca, S(k).CIhi_bca));
    fprintf('  Bootstrap BCa (FCR) CI excl. 0: %d\n', excl0(S(k).CIlo_fcr, S(k).CIhi_fcr));
    fprintf('  Bootstrap BCa (Bonf) CI excl. 0: %d\n', excl0(S(k).CIlo_bonf, S(k).CIhi_bonf));
end
fprintf('=======================================================\n');

% -------------------- Results table --------------------
Results = table(areaNames(:), ...
    meanFracU_CT1, meanFracL_CT1, log2FC_CT1, ...
    S(1).eff, S(1).CIlo_bca, S(1).CIhi_bca, S(1).CIlo_fcr, S(1).CIhi_fcr, ...
    S(1).CIlo_bonf, S(1).CIhi_bonf, S(1).perm_p, S(1).perm_q, ...
    meanFracU_CT2, meanFracL_CT2, log2FC_CT2, ...
    S(2).eff, S(2).CIlo_bca, S(2).CIhi_bca, S(2).CIlo_fcr, S(2).CIhi_fcr, ...
    S(2).CIlo_bonf, S(2).CIhi_bonf, S(2).perm_p, S(2).perm_q, ...
    S(3).eff, S(3).CIlo_bca, S(3).CIhi_bca, S(3).CIlo_fcr, S(3).CIhi_fcr, ...
    S(3).CIlo_bonf, S(3).CIhi_bonf, S(3).perm_p, S(3).perm_q, ...
    'VariableNames', {'Area', ...
    'meanFracU_CT1', 'meanFracL_CT1', 'log2FC_CT1', ...
    'Boot_logOdds_CT1', 'Boot_CIlo_BCa_CT1', 'Boot_CIhi_BCa_CT1', ...
    'Boot_CIlo_BCa_FCR_CT1', 'Boot_CIhi_BCa_FCR_CT1', ...
    'Boot_CIlo_BCa_Bonf_CT1', 'Boot_CIhi_BCa_Bonf_CT1', ...
    'Perm_p_CT1', 'Perm_q_CT1', ...
    'meanFracU_CT2', 'meanFracL_CT2', 'log2FC_CT2', ...
    'Boot_logOdds_CT2', 'Boot_CIlo_BCa_CT2', 'Boot_CIhi_BCa_CT2', ...
    'Boot_CIlo_BCa_FCR_CT2', 'Boot_CIhi_BCa_FCR_CT2', ...
    'Boot_CIlo_BCa_Bonf_CT2', 'Boot_CIhi_BCa_Bonf_CT2', ...
    'Perm_p_CT2', 'Perm_q_CT2', ...
    'Boot_Interaction', 'Boot_Int_CIlo_BCa', 'Boot_Int_CIhi_BCa', ...
    'Boot_Int_CIlo_BCa_FCR', 'Boot_Int_CIhi_BCa_FCR', ...
    'Boot_Int_CIlo_BCa_Bonf', 'Boot_Int_CIhi_BCa_Bonf', ...
    'Perm_Int_p', 'Perm_Int_q'});

% -------------------- Figures --------------------
if opt.Plot
    plot_forest_panels(Results, opt);
    plot_volcano_panels(Results, opt);
    D = struct('nu1', nu1, 'nl1', nl1, 'nu2', nu2, 'nl2', nl2, ...
        'N_U1', N_U1, 'N_L1', N_L1, 'N_U2', N_U2, 'N_L2', N_L2, ...
        'dist1', {S(1).dist}, 'dist2', {S(2).dist}, 'weighted', weighted, 'eps_prop', eps_prop);
    plot_celltype_correlation(Results, D, opt);
end
end


%% ======================= Statistics =======================
% (group_totals, group_proportion, log_odds_ratio, bca_bias_correction,
%  jackknife_acceleration, bca_ci and bh_qvalues are in functions/stats)

function S = empty_result(nAreas)
S.eff      = nan(nAreas, 1);
S.CIlo_bca = nan(nAreas, 1); S.CIhi_bca = nan(nAreas, 1);
S.CIlo     = nan(nAreas, 1); S.CIhi     = nan(nAreas, 1);
S.z0       = nan(nAreas, 1); S.a        = nan(nAreas, 1);
S.dist     = cell(nAreas, 1);
S.perm_p   = nan(nAreas, 1);
S.perm_dist = cell(nAreas, 1);
end

function S = store_result(S, i, r)
S.eff(i) = r.theta;
S.CIlo_bca(i) = r.CIlo_bca; S.CIhi_bca(i) = r.CIhi_bca;
S.CIlo(i) = r.CIlo;         S.CIhi(i) = r.CIhi;
S.z0(i) = r.z0;             S.a(i) = r.a;
S.dist{i} = r.boot;
S.perm_p(i) = r.perm_p;
S.perm_dist{i} = r.perm;
end

function p = clamp_prop(p, e)
p = max(e, min(1 - e, p));
end

function r = bias_boot_perm(sU, sL, tU, tL, nBoot, nPerm, w, e)
% Bootstrap (BCa) and permutation test of the U vs L log odds ratio
nU = numel(sU);
nL = numel(sL);

boot = zeros(nBoot, 1);
for b = 1:nBoot
    iU = randsample(nU, nU, true);
    iL = randsample(nL, nL, true);
    boot(b) = log_odds_ratio(clamp_prop(group_proportion(sU(iU), tU(iU), w), e), ...
                             clamp_prop(group_proportion(sL(iL), tL(iL), w), e));
end

pU0 = group_proportion(sU, tU, w);
pL0 = group_proportion(sL, tL, w);
theta = log_odds_ratio(pU0, pL0);
z0 = bca_bias_correction(boot, theta, nBoot);

jack = zeros(nU + nL, 1);
for j = 1:nU
    keep = setdiff(1:nU, j);
    jack(j) = log_odds_ratio(group_proportion(sU(keep), tU(keep), w), pL0);
end
for j = 1:nL
    keep = setdiff(1:nL, j);
    jack(nU + j) = log_odds_ratio(pU0, group_proportion(sL(keep), tL(keep), w));
end
a = jackknife_acceleration(jack);

[r.CIlo_bca, r.CIhi_bca] = bca_ci(boot, z0, a, 0.05);
ci = prctile(boot, [2.5, 97.5]);
r.CIlo = ci(1); r.CIhi = ci(2);

% Permutation of U/L labels
all_s = [sU; sL];
all_t = [tU; tL];
perm = zeros(nPerm, 1);
for k = 1:nPerm
    idx = randperm(nU + nL);
    perm(k) = log_odds_ratio( ...
        clamp_prop(group_proportion(all_s(idx(1:nU)), all_t(idx(1:nU)), w), e), ...
        clamp_prop(group_proportion(all_s(idx(nU+1:end)), all_t(idx(nU+1:end)), w), e));
end

r.theta = theta; r.z0 = z0; r.a = a; r.boot = boot;
r.perm = perm;
r.perm_p = sum(abs(perm) >= abs(theta)) / nPerm;
end

function r = interaction_boot_perm(s1U, s1L, s2U, s2L, t1U, t1L, t2U, t2L, nBoot, nPerm, w, e)
% Bootstrap (BCa) and permutation test of the difference in log odds ratio (CT1 - CT2)
n1U = numel(s1U); n1L = numel(s1L);
n2U = numel(s2U); n2L = numel(s2L);

boot = zeros(nBoot, 1);
for b = 1:nBoot
    i1U = randsample(n1U, n1U, true);
    i1L = randsample(n1L, n1L, true);
    i2U = randsample(n2U, n2U, true);
    i2L = randsample(n2L, n2L, true);
    lor1 = log_odds_ratio(clamp_prop(group_proportion(s1U(i1U), t1U(i1U), w), e), ...
                          clamp_prop(group_proportion(s1L(i1L), t1L(i1L), w), e));
    lor2 = log_odds_ratio(clamp_prop(group_proportion(s2U(i2U), t2U(i2U), w), e), ...
                          clamp_prop(group_proportion(s2L(i2L), t2L(i2L), w), e));
    boot(b) = lor1 - lor2;
end

p1U = group_proportion(s1U, t1U, w); p1L = group_proportion(s1L, t1L, w);
p2U = group_proportion(s2U, t2U, w); p2L = group_proportion(s2L, t2L, w);
lor1_0 = log_odds_ratio(p1U, p1L);
lor2_0 = log_odds_ratio(p2U, p2L);
theta = lor1_0 - lor2_0;
z0 = bca_bias_correction(boot, theta, nBoot);

% Jackknife: leave out one animal of any of the four groups
jack = zeros(n1U + n1L + n2U + n2L, 1);
k = 0;
for j = 1:n1U
    k = k + 1; keep = setdiff(1:n1U, j);
    jack(k) = log_odds_ratio(group_proportion(s1U(keep), t1U(keep), w), p1L) - lor2_0;
end
for j = 1:n1L
    k = k + 1; keep = setdiff(1:n1L, j);
    jack(k) = log_odds_ratio(p1U, group_proportion(s1L(keep), t1L(keep), w)) - lor2_0;
end
for j = 1:n2U
    k = k + 1; keep = setdiff(1:n2U, j);
    jack(k) = lor1_0 - log_odds_ratio(group_proportion(s2U(keep), t2U(keep), w), p2L);
end
for j = 1:n2L
    k = k + 1; keep = setdiff(1:n2L, j);
    jack(k) = lor1_0 - log_odds_ratio(p2U, group_proportion(s2L(keep), t2L(keep), w));
end
a = jackknife_acceleration(jack);

[r.CIlo_bca, r.CIhi_bca] = bca_ci(boot, z0, a, 0.05);
ci = prctile(boot, [2.5, 97.5]);
r.CIlo = ci(1); r.CIhi = ci(2);

% Permutation of cell-type labels within the U and within the L animals
all_sU = [s1U; s2U]; all_tU = [t1U; t2U];
all_sL = [s1L; s2L]; all_tL = [t1L; t2L];
perm = zeros(nPerm, 1);
for k = 1:nPerm
    iU = randperm(n1U + n2U);
    iL = randperm(n1L + n2L);
    lor1 = log_odds_ratio( ...
        clamp_prop(group_proportion(all_sU(iU(1:n1U)), all_tU(iU(1:n1U)), w), e), ...
        clamp_prop(group_proportion(all_sL(iL(1:n1L)), all_tL(iL(1:n1L)), w), e));
    lor2 = log_odds_ratio( ...
        clamp_prop(group_proportion(all_sU(iU(n1U+1:end)), all_tU(iU(n1U+1:end)), w), e), ...
        clamp_prop(group_proportion(all_sL(iL(n1L+1:end)), all_tL(iL(n1L+1:end)), w), e));
    perm(k) = lor1 - lor2;
end

r.theta = theta; r.z0 = z0; r.a = a; r.boot = boot;
r.perm = perm;
r.perm_p = sum(abs(perm) >= abs(theta)) / nPerm;
end

function ok = validateColor(x)
ok = isnumeric(x) && numel(x) == 3 && all(x >= 0) && all(x <= 1);
if ~ok, error('Color must be 1x3 RGB (0..1).'); end
end


%% ======================= Figures =======================

function plot_forest_panels(R, opt)
% Three forest plots (CT1, CT2, interaction), rows ordered by interaction effect
v1   = 100 * (0.5 * (R.meanFracU_CT1 + R.meanFracL_CT1));
v2   = 100 * (0.5 * (R.meanFracU_CT2 + R.meanFracL_CT2));
vtot = 100 * (0.25 * ((R.meanFracU_CT1 + R.meanFracL_CT1) + (R.meanFracU_CT2 + R.meanFracL_CT2)));
v1(isnan(v1)) = 0; v2(isnan(v2)) = 0; vtot(isnan(vtot)) = 0;

effTmp = R.Boot_Interaction;
effTmp(isnan(effTmp)) = -Inf;
[~, orderInt] = sort(effTmp, 'descend');

figure('Units', 'centimeters', 'Position', [2 2 34 18], 'Color', 'w');
tA = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

grey = [0.4 0.4 0.4];
panels = { ...
    R.Boot_logOdds_CT1, R.Boot_CIlo_BCa_FCR_CT1, R.Boot_CIhi_BCa_FCR_CT1, opt.ColorIntPos, v1,   'CT1: log-odds (U vs L)', 0.5; ...
    R.Boot_logOdds_CT2, R.Boot_CIlo_BCa_FCR_CT2, R.Boot_CIhi_BCa_FCR_CT2, opt.ColorIntNeg, v2,   'CT2: log-odds (U vs L)', 0.7; ...
    R.Boot_Interaction, R.Boot_Int_CIlo_BCa_FCR, R.Boot_Int_CIhi_BCa_FCR, grey,            vtot, 'CT1-CT2 interaction',    0.3};

sMin = 8; sMax = 260;
for k = 1:3
    [eff, lo, hi, col, v, ttl, textOffset] = panels{k, :};

    if max(v) == min(v)
        sVals = repmat(mean([sMin sMax]), numel(R.Area), 1);
    else
        sVals = sMin + (v - min(v)) ./ (max(v) - min(v)) * (sMax - sMin);
    end

    ax = nexttile(tA, k);
    draw_forest_aligned(ax, eff, lo, hi, R.Area, col, col, col, sVals, ttl, orderInt, opt.ForestEffectMin);
    hold(ax, 'on');

    % Marker-size legend: max, mean of positive values, min positive value
    posOnly = v(v > 0);
    if isempty(v) || all(isnan(v))
        vals = [10, 5, 1];
    elseif isempty(posOnly)
        vals = [max(v), mean(v), min(v)];
    else
        vals = [max(v), mean(posOnly), min(posOnly)];
    end
    sThree = sMin + ((vals - min(v)) / (max(v) - min(v))) * (sMax - sMin);

    yspot = [27, 28, 29];
    xlims = xlim(ax);
    for ii = 1:3
        scatter(ax, xlims(2) / 2, yspot(ii), sThree(ii), 'o', ...
            'MarkerFaceColor', opt.ColorNonSig, 'MarkerEdgeColor', opt.ColorNonSig);
        text(ax, xlims(2) / 2 + textOffset, yspot(ii), sprintf('%.1f%%', vals(ii)), ...
            'FontSize', 11, 'VerticalAlignment', 'middle');
    end
    hold(ax, 'off');
end
end

function draw_forest_aligned(ax, effects, ci_lo, ci_hi, names, colPos, colNeg, nonSigCol, svals, titleStr, orderRows, effectMin)
% Forest plot with rows in orderRows. A point is coloured colPos (colNeg)
% if effect > effectMin (< -effectMin) and its CI excludes 0.
axes(ax); hold(ax, 'on');
nTotal = numel(names);
if isempty(orderRows), orderRows = 1:nTotal; end
ypos = 1:nTotal;

allVals = [ci_lo(:); ci_hi(:); effects(:)];
finiteVals = allVals(isfinite(allVals));
if isempty(finiteVals)
    minX = -1; maxX = 1;
else
    pad = 0.08 * ((max(finiteVals) - min(finiteVals)) + eps);
    maxAbs = max(abs([min(finiteVals) - pad, max(finiteVals) + pad]));
    minX = -1.05 * maxAbs; maxX = 1.05 * maxAbs;
end

for r = 1:nTotal
    idx = orderRows(r);
    e = effects(idx); l = ci_lo(idx); h = ci_hi(idx);
    if ~isfinite(e) && ~isfinite(l) && ~isfinite(h), continue; end

    if e > effectMin && l > 0
        col = colPos;
    elseif e < -effectMin && h < 0
        col = colNeg;
    else
        col = nonSigCol;
    end

    if isfinite(l) && isfinite(h)
        plot(ax, [l h], [ypos(r) ypos(r)], '-', 'LineWidth', 1.6, 'Color', [col 0.6]);
    end
    if isfinite(e)
        scatter(ax, e, ypos(r), svals(idx), 'o', 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', col, 'MarkerFaceAlpha', 0.6);
    end
end

plot(ax, [0 0], [0 nTotal+1], '--k', 'LineWidth', 1);
set(ax, 'YTick', ypos, 'YTickLabel', names(orderRows), 'YDir', 'reverse', ...
    'FontSize', 9, 'TickDir', 'out');
xlim(ax, [minX maxX]);
xlabel(ax, 'log-odds (upper vs lower)');
title(ax, titleStr, 'FontWeight', 'bold');
box(ax, 'off');
hold(ax, 'off');
end

function plot_volcano_panels(R, opt)
% Effect vs -log10(permutation p) for CT1, CT2 and the interaction
figure('Units', 'centimeters', 'Position', [1 5 45 10], 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

panels = { ...
    R.Boot_logOdds_CT1, R.Perm_p_CT1,  R.Perm_q_CT1,  opt.ColorU1,     opt.ColorL1,     'Log-odds ratio (U vs L)',  'CT1: Volcano Plot'; ...
    R.Boot_logOdds_CT2, R.Perm_p_CT2,  R.Perm_q_CT2,  opt.ColorU2,     opt.ColorL2,     'Log-odds ratio (U vs L)',  'CT2: Volcano Plot'; ...
    R.Boot_Interaction, R.Perm_Int_p,  R.Perm_Int_q,  opt.ColorIntPos, opt.ColorIntNeg, 'Interaction (CT1 - CT2)', 'Interaction: Volcano Plot'};

for k = 1:3
    [eff, pv, qv, colPos, colNeg, xlab, ttl] = panels{k, :};
    ax = nexttile(k);
    valid = ~isnan(pv) & ~isnan(eff);
    if ~any(valid), continue; end
    effs = eff(valid); pvals = pv(valid); qvals = qv(valid); names = R.Area(valid);
    neg_log_p = -log10(pvals);

    hold on;
    for j = 1:numel(effs)
        if qvals(j) < opt.QColor(k)
            if effs(j) > 0, col = colPos; else, col = colNeg; end
        else
            col = opt.ColorNonSig;
        end
        if pvals(j) < opt.PLabel(k)
            text(effs(j) + 0.05, neg_log_p(j), names{j}, 'FontSize', 8, 'HorizontalAlignment', 'left');
        end
        scatter(effs(j), neg_log_p(j), 20, 'o', 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', col * 0.6, 'LineWidth', 0.5);
    end
    yline(-log10(0.05), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xline(0, '-k', 'LineWidth', 1);
    xlabel(xlab, 'FontSize', 10);
    ylabel('-log_{10}(p-value)', 'FontSize', 10);
    title(ttl, 'FontSize', 11);
    grid on; ax.GridAlpha = 0.1; box off;
    hold off;
end
end

function plot_celltype_correlation(R, D, opt)
% CT1 vs CT2 log odds ratios per area; error bars = MAD of the bootstrap
% distributions; r CI from a bootstrap over animals of all four groups
nAreas = numel(R.Area);
x_est = R.Boot_logOdds_CT1;
y_est = R.Boot_logOdds_CT2;
x_mad = nan(nAreas, 1);
y_mad = nan(nAreas, 1);
for i = 1:nAreas
    if ~isempty(D.dist1{i}), x_mad(i) = mad(D.dist1{i}, 1); end
    if ~isempty(D.dist2{i}), y_mad(i) = mad(D.dist2{i}, 1); end
end

validPair = isfinite(x_est) & isfinite(y_est);
if sum(validPair) <= 2, return; end
xv = x_est(validPair); yv = y_est(validPair);
xm = x_mad(validPair); ym = y_mad(validPair);
[r_obs, p_obs] = corr(xv, yv, 'type', 'Pearson');

nU1 = size(D.nu1, 1); nL1 = size(D.nl1, 1);
nU2 = size(D.nu2, 1); nL2 = size(D.nl2, 1);
nboot = opt.NBootstrap;
boot_r = zeros(nboot, 1);
x_fit = linspace(min(xv), max(xv), 100);
y_boot = zeros(nboot, numel(x_fit));

rng(opt.Seed);
for b = 1:nboot
    iU1 = randi(nU1, nU1, 1);
    iL1 = randi(nL1, nL1, 1);
    iU2 = randi(nU2, nU2, 1);
    iL2 = randi(nL2, nL2, 1);
    lor1 = nan(nAreas, 1);
    lor2 = nan(nAreas, 1);
    for i = 1:nAreas
        lor1(i) = resampled_log_odds(D.nu1(iU1, i), D.nl1(iL1, i), D.N_U1(iU1), D.N_L1(iL1), D);
        lor2(i) = resampled_log_odds(D.nu2(iU2, i), D.nl2(iL2, i), D.N_U2(iU2), D.N_L2(iL2), D);
    end
    ok = isfinite(lor1) & isfinite(lor2);
    if sum(ok) > 2
        boot_r(b) = corr(lor1(ok), lor2(ok), 'type', 'Pearson');
        y_boot(b, :) = polyval(polyfit(lor1(ok), lor2(ok), 1), x_fit);
    else
        boot_r(b) = NaN;
        y_boot(b, :) = NaN;
    end
end

ci_r = prctile(boot_r, [2.5 97.5]);
y_lo = prctile(y_boot, 2.5, 1);
y_hi = prctile(y_boot, 97.5, 1);
y_med = median(y_boot, 1, 'omitnan');

figure('Position', [329.8 326.6 482.4 344.8]);
hold on;
fill([x_fit, fliplr(x_fit)], [y_lo, fliplr(y_hi)], [0.3 0.3 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(x_fit, y_med, 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
errorbar(xv, yv, ym, ym, xm, xm, 'Color', [0.4 0.4 0.4], 'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 0.9);
scatter(xv, yv, 30, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 1);
set(gca, 'FontSize', 12);
xlabel('log-odds (upper/lower) CT1');
ylabel('log-odds (upper/lower) CT2');
tit = title(sprintf('r=%.2f [%.2f, %.2f], p=%.3g', r_obs, ci_r(1), ci_r(2), p_obs));
tit.FontWeight = 'normal';
tit.FontSize = 12;

figure;
histogram(boot_r, 50, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.6]);
hold on;
xline(r_obs, 'r-', 'LineWidth', 2, 'Label', sprintf('r=%.2f', r_obs));
xline(ci_r(1), 'b--', 'LineWidth', 1.5);
xline(ci_r(2), 'b--', 'LineWidth', 1.5);
xlabel('Pearson r');
ylabel('Probability');
title(sprintf('Bootstrap Distribution: r=%.2f [%.2f, %.2f]', r_obs, ci_r(1), ci_r(2)));
set(gca, 'FontSize', 12);

fprintf('\nLog-odds correlation (CT1 vs CT2):\n');
fprintf('  Observed r: %.3f\n', r_obs);
fprintf('  Bootstrap mean r: %.3f\n', mean(boot_r, 'omitnan'));
fprintf('  Bootstrap median r: %.3f\n', median(boot_r, 'omitnan'));
fprintf('  Bootstrap std: %.3f\n', std(boot_r, 'omitnan'));
end

function lor = resampled_log_odds(sU, sL, tU, tL, D)
% Log odds ratio of one resample; NaN if the area has no data
lor = NaN;
if any(isnan(sU)) || any(isnan(sL)), return; end
pU = clamp_prop(group_proportion(sU, tU, D.weighted), D.eps_prop);
pL = clamp_prop(group_proportion(sL, tL, D.weighted), D.eps_prop);
lor = log(pU / (1 - pU)) - log(pL / (1 - pL));
end
