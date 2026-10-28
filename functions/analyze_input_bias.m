function Results = analyze_input_bias(nu, nl, areaNames, varargin)
%ANALYZE_INPUT_BIAS Test each area for an input bias between two injection groups.
%   Results = ANALYZE_INPUT_BIAS(nu, nl, areaNames, Name, Value, ...)
%
%   nu, nl     animals-by-areas cell counts for the upper (U, medial sSC)
%              and lower (L, lateral sSC) visual field injections
%   areaNames  area acronyms, one per column of nu and nl
%
%   For each area the effect is the log odds ratio of the input fraction,
%   U vs L (> 0: bias towards the upper visual field). Input fractions are
%   counts / TotalU (TotalL) per animal. The effect is tested with
%     - a bootstrap over animals (BCa CIs: unadjusted, Bonferroni, FCR)
%     - a permutation test of the group labels (p-values; CIs by test inversion)
%     - a binomial GLM, with an observation-level random effect added when
%       the data are overdispersed (Pearson dispersion > 1.5)
%     - a rank-sum test and a t-test on the fractions
%
%   Name-value options [default]
%     TotalU, TotalL        total inputs per animal [sum(nu,2), sum(nl,2)]
%     Weighted              false: mean of per-animal proportions;
%                           true: pooled proportions [false]
%     CorrectContamination  un-mix proportions using the injection spread [false]
%     EpsU, EpsL            fraction of the U (L) injection spreading into
%                           L (U) territory [0.033, 0.162]
%     NBootstrap            bootstrap resamples [10000]
%     NPermutations         permutations [10000]
%     Seed                  random seed [42]
%     Alpha                 pseudocount for log2FC and logit [1e-10]
%     QDiscoveryPerm        q cutoff counting discoveries for the permutation
%                           FCR-adjusted CIs [0.045]
%     QDiscoveryBoot        q cutoff counting discoveries for the bootstrap
%                           FCR-adjusted CIs [0.035]
%     QPlot                 q cutoff for colouring areas as significant [0.046]
%     Plot                  draw figures [true]
%     Figures               figures to draw, "all" or a list of: "bootForest",
%                           "permVolcano", "permPvalues", "permBootCI",
%                           "glmmVolcano", "glmmForestVertical",
%                           "glmmForest", "glmmBoxplots", "permInvertedCI" ["all"]
%     ColorU, ColorL, ColorNonSig, PosPalette, NegPalette   plot colours
%     plot_seq              order of areas in the plots [sorted by effect]
%     GroupNames            names of the two groups in axis labels ["upper", "lower"]
%
%   Results has one row per area. Boot_q holds the Benjamini-Hochberg
%   q-values of the PERMUTATION p-values (used for the bootstrap FCR CIs).
%   GLMM_CIlo/GLMM_CIhi are Bonferroni-adjusted (same as *_Bonf).

% -------------------- Options --------------------
p = inputParser;
addRequired(p, 'nu', @isnumeric);
addRequired(p, 'nl', @isnumeric);
addOptional(p, 'areaNames', [], @(x) iscellstr(x) || isstring(x) || isempty(x));
addParameter(p, 'TotalU', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'TotalL', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'Weighted', false, @islogical);
addParameter(p, 'CorrectContamination', false, @islogical);
addParameter(p, 'EpsU', 0.033, @isnumeric);
addParameter(p, 'EpsL', 0.162, @isnumeric);
addParameter(p, 'NBootstrap', 10000, @isnumeric);
addParameter(p, 'NPermutations', 10000, @isnumeric);
addParameter(p, 'Seed', 42, @isnumeric);
addParameter(p, 'Alpha', 1e-10, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'QDiscoveryPerm', 0.045, @isnumeric);
addParameter(p, 'QDiscoveryBoot', 0.035, @isnumeric);
addParameter(p, 'QPlot', 0.046, @isnumeric);
addParameter(p, 'Plot', true, @islogical);
addParameter(p, 'Figures', "all", @(x) isstring(x) || iscellstr(x) || ischar(x));
addParameter(p, 'ColorNonSig', [0.7 0.7 0.7], @validateColor);
addParameter(p, 'PosPalette', [0.8510 0.5000 0; 0.8808 0.6000 0.2000; 0.9106 0.7000 0.4000; ...
    0.9404 0.8000 0.6000; 0.9702 0.9000 0.8000], @validatePalette);
addParameter(p, 'NegPalette', [0.2100 0.3900 0.7300; 0.3680 0.5120 0.7840; 0.5260 0.6340 0.8380; ...
    0.6840 0.7560 0.8920; 0.8420 0.8780 0.9460], @validatePalette);
addParameter(p, 'ColorU', [0.2 0.4 0.8], @validateColor);
addParameter(p, 'ColorL', [0.85 0.2 0.2], @validateColor);
addParameter(p, 'plot_seq', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'GroupNames', ["upper", "lower"], @(x) (isstring(x) || iscellstr(x)) && numel(x) == 2);
parse(p, nu, nl, areaNames, varargin{:});
opt = p.Results;

weighted   = opt.Weighted;
correct    = opt.CorrectContamination;
epsU       = opt.EpsU;
epsL       = opt.EpsL;
nBootstrap = opt.NBootstrap;
nPerm      = opt.NPermutations;
alpha      = opt.Alpha;
eps_prop   = 1e-10;            % keeps proportions inside (0, 1)

% -------------------- Sizes, names, totals --------------------
[nU_anim, nAreas] = size(nu);
[nL_anim, nAreasL] = size(nl);
assert(nAreasL == nAreas, 'nu and nl must have the same number of columns (areas).');
if isempty(opt.areaNames)
    areaNames = arrayfun(@(k) sprintf('Area_%02d', k), 1:nAreas, 'UniformOutput', false);
else
    areaNames = cellstr(opt.areaNames);
end

N_U = group_totals(opt.TotalU, nu, nU_anim, 'TotalU');
N_L = group_totals(opt.TotalL, nl, nL_anim, 'TotalL');
if any(N_U == 0) || any(N_L == 0), error('TotalU/TotalL must be > 0'); end

% -------------------- Fractions --------------------
fracU = nu ./ N_U;               % per-animal input fractions
fracL = nl ./ N_L;
meanFracU = mean(fracU, 1);
meanFracL = mean(fracL, 1);
log2FC = log2((meanFracU + alpha) ./ (meanFracL + alpha));

if correct
    det_m = (1 - epsU) * (1 - epsL) - epsU * epsL;
    tu = ((1 - epsL) * meanFracU - epsU * meanFracL) / det_m;
    tl = ((1 - epsU) * meanFracL - epsL * meanFracU) / det_m;
    fprintf('Contamination correction ON (eps_U=%.3f, eps_L=%.3f, det=%.3f).\n', epsU, epsL, det_m);
    fprintf('  Areas clamped after un-mixing: %d of %d\n', sum(tu < 0 | tl < 0), numel(tu));
end

% -------------------- Storage --------------------
ranksum_p = nan(nAreas, 1);  ttest_p = nan(nAreas, 1);
GLMM_effect = nan(nAreas, 1); GLMM_SE = nan(nAreas, 1); GLMM_p = nan(nAreas, 1);

bootstrap_effect = nan(nAreas, 1);
bootstrap_CIlo = nan(nAreas, 1);      bootstrap_CIhi = nan(nAreas, 1);
bootstrap_CIlo_bca = nan(nAreas, 1);  bootstrap_CIhi_bca = nan(nAreas, 1);
bootstrap_p = nan(nAreas, 1);
bootstrap_z0 = nan(nAreas, 1);        bootstrap_a = nan(nAreas, 1);
bootstrap_distributions = cell(nAreas, 1);

permutation_effect = nan(nAreas, 1);
permutation_p = nan(nAreas, 1);
permutation_CIlo = nan(nAreas, 1);    permutation_CIhi = nan(nAreas, 1);
permutation_distributions = cell(nAreas, 1);

% -------------------- Per-area tests --------------------
% Random numbers are drawn in a fixed order (bootstrap, then permutations,
% area by area), so results are reproducible for a given seed.
rng(opt.Seed);

trials_U = N_U(:);
trials_L = N_L(:);

for i = 1:nAreas
    dataU = fracU(:, i);
    dataL = fracL(:, i);
    succ_U = nu(:, i);
    succ_L = nl(:, i);
    nU = numel(succ_U);
    nL = numel(succ_L);

    % --- Binomial GLM(M) ---
    [GLMM_effect(i), GLMM_SE(i), GLMM_p(i)] = fit_binomial_glmm(succ_U, succ_L, trials_U, trials_L, i);

    % --- BCa bootstrap over animals ---
    try
        boot = zeros(nBootstrap, 1);
        for b = 1:nBootstrap
            iU = randsample(nU, nU, true);
            iL = randsample(nL, nL, true);
            pU = group_proportion(succ_U(iU), trials_U(iU), weighted);
            pL = group_proportion(succ_L(iL), trials_L(iL), weighted);
            pU = max(eps_prop, min(1 - eps_prop, pU));
            pL = max(eps_prop, min(1 - eps_prop, pL));
            [pU, pL, ok] = unmix(pU, pL, epsU, epsL, correct);
            if ~ok
                boot(b) = NaN;
                continue;
            end
            boot(b) = log_odds_ratio(pU, pL);
        end

        nDropped = sum(isnan(boot));
        if nDropped > 0
            fprintf('  Area %d: %d/%d bootstrap resamples dropped\n', i, nDropped, nBootstrap);
        end
        boot = boot(~isnan(boot));

        % Observed effect
        [pU0, pL0] = unmix(group_proportion(succ_U, trials_U, weighted), ...
                           group_proportion(succ_L, trials_L, weighted), epsU, epsL, correct);
        theta_hat = log_odds_ratio(pU0, pL0);

        z0 = bca_bias_correction(boot, theta_hat, nBootstrap);

        % Acceleration (jackknife, leaving out one animal of either group)
        jack = zeros(nU + nL, 1);
        for j = 1:nU
            keep = setdiff(1:nU, j);
            [a1, b1] = unmix(group_proportion(succ_U(keep), trials_U(keep), weighted), ...
                             group_proportion(succ_L, trials_L, weighted), epsU, epsL, correct);
            jack(j) = log_odds_ratio(a1, b1);
        end
        for j = 1:nL
            keep = setdiff(1:nL, j);
            [a1, b1] = unmix(group_proportion(succ_U, trials_U, weighted), ...
                             group_proportion(succ_L(keep), trials_L(keep), weighted), epsU, epsL, correct);
            jack(nU + j) = log_odds_ratio(a1, b1);
        end
        a = jackknife_acceleration(jack);

        % Confidence intervals (unadjusted)
        [bootstrap_CIlo_bca(i), bootstrap_CIhi_bca(i)] = bca_ci(boot, z0, a, 0.05);
        ci = prctile(boot, [2.5, 97.5]);
        bootstrap_CIlo(i) = ci(1);
        bootstrap_CIhi(i) = ci(2);

        bootstrap_effect(i) = theta_hat;
        bootstrap_z0(i) = z0;
        bootstrap_a(i) = a;

        % Two-sided bootstrap p-value
        pb = sum(boot > 0) / nBootstrap;
        if pb > 0.5, pb = 1 - pb; end
        bootstrap_p(i) = 2 * pb;

        bootstrap_distributions{i} = boot;
    catch ME
        warning('Bootstrap failed for area %d: %s', i, ME.message);
        bootstrap_effect(i) = NaN;
        bootstrap_CIlo(i) = NaN;     bootstrap_CIhi(i) = NaN;
        bootstrap_CIlo_bca(i) = NaN; bootstrap_CIhi_bca(i) = NaN;
        bootstrap_p(i) = NaN;
        bootstrap_z0(i) = NaN;       bootstrap_a(i) = NaN;
        bootstrap_distributions{i} = [];
    end

    % --- Permutation test of group labels ---
    try
        all_succ   = [succ_U; succ_L];
        all_trials = [trials_U; trials_L];
        nTot = nU + nL;

        [pU0, pL0] = unmix(group_proportion(succ_U, trials_U, weighted), ...
                           group_proportion(succ_L, trials_L, weighted), epsU, epsL, correct);
        pU0 = max(eps_prop, min(1 - eps_prop, pU0));
        pL0 = max(eps_prop, min(1 - eps_prop, pL0));
        theta_hat = log_odds_ratio(pU0, pL0);
        permutation_effect(i) = theta_hat;

        perm = zeros(nPerm, 1);
        for k = 1:nPerm
            idx = randperm(nTot);
            pU = group_proportion(all_succ(idx(1:nU)), all_trials(idx(1:nU)), weighted);
            pL = group_proportion(all_succ(idx(nU+1:end)), all_trials(idx(nU+1:end)), weighted);
            pU = max(eps_prop, min(1 - eps_prop, pU));
            pL = max(eps_prop, min(1 - eps_prop, pL));
            [pU, pL] = unmix(pU, pL, epsU, epsL, correct);
            perm(k) = log_odds_ratio(pU, pL);
        end

        permutation_p(i) = sum(abs(perm) >= abs(theta_hat)) / nPerm;
        [permutation_CIlo(i), permutation_CIhi(i)] = inverted_perm_ci(perm, theta_hat, 0.05, nPerm);
        permutation_distributions{i} = perm;
    catch ME
        warning('Permutation test failed for area %d: %s', i, ME.message);
        permutation_effect(i) = NaN;
        permutation_p(i) = NaN;
        permutation_CIlo(i) = NaN;
        permutation_CIhi(i) = NaN;
        permutation_distributions{i} = [];
    end

    % --- Rank-sum and t-test on the fractions ---
    try ranksum_p(i) = ranksum(dataU, dataL); catch, ranksum_p(i) = NaN; end
    try [~, ttest_p(i)] = ttest2(dataU, dataL); catch, ttest_p(i) = NaN; end
end

% -------------------- Multiple-comparison adjustments --------------------
qvals_perm = bh_qvalues(permutation_p);
ranksum_q  = bh_qvalues(ranksum_p);
qvals      = bh_qvalues(GLMM_p);                       % GLMM q-values

% Permutation CIs: FCR-adjusted (discoveries at q < QDiscoveryPerm) and Bonferroni
validIdx_perm = ~isnan(permutation_p);
m_perm = sum(validIdx_perm);
R_perm = sum(qvals_perm(validIdx_perm) < opt.QDiscoveryPerm);
if R_perm == 0
    alpha_adj_perm = 0.05 / m_perm;
    fprintf('Permutation FCR: No discoveries. Using Bonferroni (alpha = %.5f).\n', alpha_adj_perm);
else
    alpha_adj_perm = 0.05 * (R_perm / m_perm);
    fprintf('Permutation FCR: %d discoveries. Adjusted alpha = %.5f.\n', R_perm, alpha_adj_perm);
end

permutation_CIlo_fcr  = nan(nAreas, 1); permutation_CIhi_fcr  = nan(nAreas, 1);
permutation_CIlo_bonf = nan(nAreas, 1); permutation_CIhi_bonf = nan(nAreas, 1);
for i = 1:nAreas
    if ~isempty(permutation_distributions{i}) && ~isnan(permutation_effect(i))
        [permutation_CIlo_fcr(i), permutation_CIhi_fcr(i)] = inverted_perm_ci( ...
            permutation_distributions{i}, permutation_effect(i), alpha_adj_perm, nPerm);
        [permutation_CIlo_bonf(i), permutation_CIhi_bonf(i)] = inverted_perm_ci( ...
            permutation_distributions{i}, permutation_effect(i), 0.05 / nAreas, nPerm);
    end
end

% GLMM CIs: Bonferroni
z_bonf = -norminv((0.05 / nAreas) / 2);
GLMM_CIlo_bonf = GLMM_effect - z_bonf * GLMM_SE;
GLMM_CIhi_bonf = GLMM_effect + z_bonf * GLMM_SE;

% Bootstrap CIs: FCR-adjusted (discoveries counted on the PERMUTATION
% q-values at q < QDiscoveryBoot) and Bonferroni
m_boot = sum(~isnan(bootstrap_p));
qvals_boot = nan(nAreas, 1);
if m_boot > 0
    qvals_boot = qvals_perm;
end
R_boot = sum(qvals_boot(validIdx_perm) < opt.QDiscoveryBoot);
if R_boot == 0
    alpha_adj_boot = 0.05 / m_boot;
    fprintf('Bootstrap FCR: No discoveries. Using Bonferroni.\n');
else
    alpha_adj_boot = 0.05 * (R_boot / m_boot);
    fprintf('Bootstrap FCR: %d discoveries. Adjusting Alpha to %.5f.\n', R_boot, alpha_adj_boot);
end

bootstrap_CIlo_bca_fcr  = nan(nAreas, 1); bootstrap_CIhi_bca_fcr  = nan(nAreas, 1);
bootstrap_CIlo_bca_bonf = nan(nAreas, 1); bootstrap_CIhi_bca_bonf = nan(nAreas, 1);
for i = 1:nAreas
    if ~isempty(bootstrap_distributions{i}) && ~isnan(bootstrap_z0(i))
        [bootstrap_CIlo_bca_fcr(i), bootstrap_CIhi_bca_fcr(i)] = bca_ci( ...
            bootstrap_distributions{i}, bootstrap_z0(i), bootstrap_a(i), alpha_adj_boot);
        [bootstrap_CIlo_bca_bonf(i), bootstrap_CIhi_bca_bonf(i)] = bca_ci( ...
            bootstrap_distributions{i}, bootstrap_z0(i), bootstrap_a(i), 0.05 / nAreas);
    end
end

fprintf('Permutation: Calculated FDR-adjusted q-values for %d tests.\n', m_perm);

% -------------------- Results table --------------------
Results = table(areaNames(:), meanFracU(:), meanFracL(:), log2FC(:), ...
    GLMM_effect, GLMM_SE, GLMM_CIlo_bonf, GLMM_CIhi_bonf, ...
    GLMM_CIlo_bonf, GLMM_CIhi_bonf, ...
    GLMM_p, qvals, ...
    bootstrap_effect, bootstrap_CIlo, bootstrap_CIhi, ...
    bootstrap_CIlo_bca, bootstrap_CIhi_bca, ...
    bootstrap_CIlo_bca_bonf, bootstrap_CIhi_bca_bonf, ...
    bootstrap_CIlo_bca_fcr, bootstrap_CIhi_bca_fcr, ...
    bootstrap_p, qvals_boot, bootstrap_z0, bootstrap_a, ...
    permutation_effect, permutation_p, qvals_perm, permutation_CIlo, permutation_CIhi, ...
    permutation_CIlo_bonf, permutation_CIhi_bonf, permutation_CIlo_fcr, permutation_CIhi_fcr, ...
    ranksum_p, ranksum_q, ttest_p, ...
    'VariableNames', {'Area', 'meanFracU', 'meanFracL', 'log2FC', ...
    'GLMM_logOdds_U_vs_L', 'GLMM_SE', 'GLMM_CIlo', 'GLMM_CIhi', ...
    'GLMM_CIlo_Bonf', 'GLMM_CIhi_Bonf', ...
    'GLMM_p', 'GLMM_q', ...
    'Boot_logOdds_U_vs_L', 'Boot_CIlo_percentile', 'Boot_CIhi_percentile', ...
    'Boot_CIlo_BCa', 'Boot_CIhi_BCa', ...
    'Boot_CIlo_BCa_Bonf', 'Boot_CIhi_BCa_Bonf', ...
    'Boot_CIlo_BCa_FCR', 'Boot_CIhi_BCa_FCR', ...
    'Boot_p', 'Boot_q', 'Boot_z0', 'Boot_a', ...
    'Perm_logOdds_U_vs_L', 'Perm_p', 'Perm_q', 'Perm_CIlo', 'Perm_CIhi', ...
    'Perm_CIlo_Bonf', 'Perm_CIhi_Bonf', 'Perm_CIlo_FCR', 'Perm_CIhi_FCR', ...
    'ranksum_p', 'ranksum_q', 'ttest_p'});

Results.Bootstrap_Distributions   = bootstrap_distributions;
Results.Permutation_Distributions = permutation_distributions;

% -------------------- Figures --------------------
if opt.Plot
    pl = opt;
    pl.nAreas = nAreas;
    pl.fracU = fracU;
    pl.fracL = fracL;
    pl.GroupNames = cellstr(opt.GroupNames);
    if isequal(pl.GroupNames, {'upper', 'lower'})
        pl.GroupTicks = {'U', 'L'};
    else
        pl.GroupTicks = pl.GroupNames;
    end
    figs = string(opt.Figures);
    draw = @(name) any(figs == "all") || any(figs == name);

    if draw("bootForest"),         plot_boot_forest(Results, pl); end
    if draw("permVolcano"),        plot_perm_volcano(Results, pl); end
    if draw("permPvalues"),        plot_perm_pvalues(Results, pl); end
    if draw("permBootCI"),         plot_perm_boot_ci(Results, pl); end
    if draw("glmmVolcano"),        plot_glmm_volcano(Results, pl); end
    if draw("glmmForestVertical"), plot_glmm_forest_vertical(Results, pl); end
    if draw("glmmForest"),         plot_glmm_forest(Results, pl); end
    if draw("glmmBoxplots"),       plot_glmm_boxplots(Results, pl); end
    if draw("permInvertedCI"),     plot_perm_inverted_ci(Results, pl); end
end

fprintf('Notes: q = Benjamini-Hochberg FDR-adjusted p-value.\n');
end

%% ======================= Statistics helpers =======================
% (group_totals, group_proportion, log_odds_ratio, bca_bias_correction,
%  jackknife_acceleration, bca_ci and bh_qvalues are in functions/stats)

function [lo, hi] = inverted_perm_ci(permDist, theta_hat, alphaLevel, nPerm)
% Shifts that are not rejected by the permutation test at alphaLevel
delta_grid = linspace(-6, 6, 150);
p_shifted = sum(abs(permDist) >= abs(theta_hat - delta_grid), 1) / nPerm;
not_rejected = delta_grid(~(p_shifted < alphaLevel));
if isempty(not_rejected)
    lo = theta_hat;
    hi = theta_hat;
else
    lo = min(not_rejected);
    hi = max(not_rejected);
end
end

function [effect, se, pval] = fit_binomial_glmm(succ_U, succ_L, trials_U, trials_L, areaIdx)
% Binomial GLM of U vs L; observation-level random effect if overdispersed
effect = NaN; se = NaN; pval = NaN;
nU = numel(succ_U);
nL = numel(succ_L);
Successes = [succ_U; succ_L];
Trials    = [trials_U; trials_L];
Group     = categorical([repmat("U", nU, 1); repmat("L", nL, 1)], ["U", "L"]);
tbl = table(Successes, Trials, Group);
try
    glme = fitglme(tbl, 'Successes ~ Group', 'Distribution', 'Binomial', ...
        'BinomialSize', tbl.Trials, 'FitMethod', 'Laplace', 'Verbose', 0);

    % Pearson dispersion
    p_hat   = predict(glme);
    resid   = Successes - p_hat .* Trials;
    var_hat = Trials .* p_hat .* (1 - p_hat);
    nObs = numel(Successes);
    nFix = height(glme.Coefficients);
    dispersion = sum(resid.^2 ./ max(var_hat, eps)) / max(nObs - nFix, 1);

    if dispersion > 1.5
        tbl.Obs = categorical(string((1:nObs)'));
        glme = fitglme(tbl, 'Successes ~ Group + (1|Obs)', 'Distribution', 'Binomial', ...
            'BinomialSize', tbl.Trials, 'FitMethod', 'Laplace', 'Verbose', 0);
    end

    coefTbl = glme.Coefficients;
    idx = strcmp(coefTbl.Name, 'Group_L');
    if ~any(idx), idx = contains(coefTbl.Name, 'Group'); end
    assert(nnz(idx) == 1, 'Group coefficient not found.');
    effect = -coefTbl.Estimate(idx);     % sign flipped: U vs L
    se     = coefTbl.SE(idx);
    pval   = coefTbl.pValue(idx);
catch ME
    warning('GLME failed for area %d: %s', areaIdx, ME.message);
    effect = NaN; se = NaN; pval = NaN;
end
end

function [pU, pL, ok] = unmix(pU, pL, eU, eL, doIt)
% Correct U and L proportions for injection spread into the other territory
ok = true;
if ~doIt, return; end
det_m = (1 - eU) * (1 - eL) - eU * eL;
if abs(det_m) < 1e-6, return; end
pU_new = ((1 - eL) * pU - eU * pL) / det_m;
pL_new = ((1 - eU) * pL - eL * pU) / det_m;
if pU_new <= 0 || pL_new <= 0
    ok = false;
    return;
end
pU = pU_new;
pL = pL_new;
end

function ok = validateColor(x)
ok = isnumeric(x) && isequal(size(x), [1, 3]) && all(x >= 0) && all(x <= 1);
if ~ok, error('Color inputs must be 1x3 RGB in 0..1'); end
end

function ok = validatePalette(x)
ok = isnumeric(x) && size(x, 2) == 3 && all(x(:) >= 0) && all(x(:) <= 1);
if ~ok, error('Palette must be an Mx3 numeric matrix with values in 0..1'); end
end

%% ======================= Plot helpers =======================

function [vals, sidx] = sort_for_plot(vals, plot_seq, direction)
% Order given by plot_seq, otherwise sorted by value
if isempty(plot_seq)
    [vals, sidx] = sort(vals, direction);
else
    sidx = plot_seq;
    vals = vals(plot_seq);
end
end

function [sVals, legVals, legSizes] = abundance_marker_sizes(abund_pct)
% Marker areas scaled linearly with input % (plus a small offset)
sMin = 5; sMax = 225;
posVals = abund_pct(abund_pct > 0);
if isempty(posVals)
    offset = 0.1;
else
    offset = max(min(posVals) / 10, 0.01);
end
Vals = abund_pct + offset;
minV = min(Vals);
denom = max(Vals) - minV;
if denom <= eps
    sVals = repmat(mean([sMin sMax]), numel(abund_pct), 1);
else
    sVals = sMin + ((Vals - minV) / denom) * (sMax - sMin);
end

% Legend: max, mean of positive values, min positive value
if isempty(abund_pct) || all(isnan(abund_pct))
    legVals = [10, 5, 1];
else
    posOnly = abund_pct(abund_pct > 0);
    if isempty(posOnly)
        legVals = [max(abund_pct), mean(abund_pct), min(abund_pct)];
    else
        legVals = [max(abund_pct), mean(posOnly), min(posOnly)];
    end
end
if denom <= eps
    legSizes = repmat(mean([sMin sMax]), 1, 3);
else
    legSizes = sMin + (((legVals + offset) - minV) / denom) * (sMax - sMin);
end
end

function col = palette_color(mag, palette)
% Colour from palette; mag = 1 gives the first (darkest) entry
m = size(palette, 1);
xi = linspace(0, 1, m);
col = interp1(xi, palette, 1 - mag, 'linear');
col = min(max(col, 0), 1);
end

function cols = significance_colors(effs, qvals, pl)
% ColorU / ColorL for q < QPlot (by sign of effect), grey otherwise
neutralDot = [0.7 0.7 0.7];
cols = repmat(neutralDot, numel(effs), 1);
for k = 1:numel(effs)
    if effs(k) > 0 && qvals(k) < pl.QPlot
        cols(k, :) = pl.ColorU;
    elseif effs(k) <= 0 && qvals(k) < pl.QPlot
        cols(k, :) = pl.ColorL;
    end
end
end

function draw_forest_columns(eff, lo, hi, names, sVals, cols, lineCols, xMargin, ...
    legVals, legSizes, legStep, ylab, ttl)
% Effects with CIs, one column per area; marker size = input %
neutralDot = [0.7 0.7 0.7];
nplot = numel(eff);
xpos = 1:nplot;
hold on;
for k = 1:nplot
    plot([xpos(k) xpos(k)], [lo(k) hi(k)], '-', 'LineWidth', 1.8, 'Color', lineCols(k, :));
    scatter(xpos(k), eff(k), sVals(k), 'o', 'MarkerFaceColor', cols(k, :), ...
        'MarkerEdgeColor', cols(k, :) * 0.5, 'LineWidth', 0.5, 'MarkerFaceAlpha', 1);
end
plot([0 nplot+1], [0 0], '--k', 'LineWidth', 1);

pad = 0.5;
ylim([min(lo) - pad, max(hi) + pad]);
xlim([1 - xMargin, nplot + xMargin]);
set(gca, 'XTick', xpos, 'XTickLabel', names, 'FontSize', 10, ...
    'TickDir', 'out', 'XTickLabelRotation', 90);
ax = gca;
box off;
set(gca, 'YGrid', 'on', 'XGrid', 'off');
ax.GridColor = [0.2 0.2 0.2];
ax.GridAlpha = 0.05;
ylabel(ylab, 'FontSize', 11);
if ~isempty(ttl), title(ttl, 'FontSize', 12); end

% Marker-size legend
ylims = ylim;
yspot = ylims(2) - 0.2 - [0, legStep, 2 * legStep];
xlims = xlim;
for i = 1:3
    scatter(gca, xlims(2) - 3, yspot(i), legSizes(i), 'o', ...
        'MarkerFaceColor', neutralDot, 'MarkerEdgeColor', neutralDot);
    text(gca, xlims(2) - 2.5, yspot(i), sprintf('%.1f%%', legVals(i)), ...
        'FontSize', 11, 'VerticalAlignment', 'middle');
end
hold off;
end

%% ======================= Figures =======================

function plot_boot_forest(R, pl)
% Bootstrap effects with FCR-adjusted BCa CIs, coloured by effect size
figure('Units', 'centimeters', 'Position', [15.5998 10.7103 21.2513 5.7362], 'Color', 'w');
valid = ~isnan(R.Boot_logOdds_U_vs_L);
if ~any(valid)
    text(0.1, 0.5, 'No valid Bootstrap results to plot', 'FontSize', 11);
    return;
end
effs  = R.Boot_logOdds_U_vs_L(valid);
ci_lo = R.Boot_CIlo_BCa_FCR(valid);
ci_hi = R.Boot_CIhi_BCa_FCR(valid);
names = R.Area(valid);
abund_pct = 100 * (R.meanFracU(valid) + R.meanFracL(valid)) ./ 2;

[eff_sorted, sidx] = sort_for_plot(effs, pl.plot_seq, 'descend');
ci_lo = ci_lo(sidx); ci_hi = ci_hi(sidx); names = names(sidx); abund_pct = abund_pct(sidx);

[sVals, legVals, legSizes] = abundance_marker_sizes(abund_pct);

% Colour scaled within the range of positive (negative) effects
pos_effects = eff_sorted(eff_sorted > 0);
neg_effects = eff_sorted(eff_sorted <= 0);
if isempty(pos_effects), maxPos = eps; minPos = 0;
else, maxPos = max(pos_effects); minPos = min(pos_effects); end
if isempty(neg_effects), maxNeg = eps; minNeg = 0;
else, maxNeg = max(abs(neg_effects)); minNeg = min(abs(neg_effects)); end

nplot = numel(eff_sorted);
cols = repmat([0.7 0.7 0.7], nplot, 1);
for k = 1:nplot
    e = eff_sorted(k);
    if e > 0
        if maxPos > minPos, mag = (e - minPos) / (maxPos - minPos); else, mag = 0.5; end
        cols(k, :) = palette_color(max(0, min(1, mag)), pl.PosPalette);
    elseif e <= 0
        if maxNeg > minNeg, mag = (abs(e) - minNeg) / (maxNeg - minNeg); else, mag = 0.5; end
        cols(k, :) = palette_color(max(0, min(1, mag)), pl.NegPalette);
    end
end
lineCols = max(0, cols * 0.9);

draw_forest_columns(eff_sorted, ci_lo, ci_hi, names, sVals, cols, lineCols, 0.5, ...
    legVals, legSizes, 0.6, sprintf('btsrp log-odds (%s vs %s)', pl.GroupNames{:}), '');
end

function plot_perm_volcano(R, pl)
% Effect vs -log10(permutation p); coloured if q < QPlot, labelled if p < 0.05
figure('Units', 'centimeters', 'Position', [6.7098 10.3505 10.3505 8.0010], 'Color', 'w');
valid = ~isnan(R.Perm_p) & ~isnan(R.Perm_logOdds_U_vs_L);
if ~any(valid), return; end
effs  = R.Perm_logOdds_U_vs_L(valid);
pvals = R.Perm_p(valid);
qvals = R.Perm_q(valid);
names = R.Area(valid);
neg_log_p = -log10(pvals);

hold on;
for k = 1:numel(effs)
    if qvals(k) < pl.QPlot
        if effs(k) > 0, col = pl.ColorU; else, col = pl.ColorL; end
    else
        col = pl.ColorNonSig;
    end
    scatter(effs(k), neg_log_p(k), 20, 'o', 'MarkerFaceColor', col, ...
        'MarkerEdgeColor', col * 0.6, 'LineWidth', 0.5);
    if pvals(k) < 0.05
        text(effs(k) + 0.05, neg_log_p(k), names{k}, 'FontSize', 9, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end
yline(-log10(0.05), '--', 'p = 0.05', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xline(0, '-k', 'LineWidth', 1);
xlabel(sprintf('Log-odds ratio (%s vs %s)', pl.GroupNames{:}), 'FontSize', 11);
ylabel('-log_{10}(p-value)', 'FontSize', 11);
title('Permutation Test: Volcano Plot', 'FontSize', 12);
grid on;
ax = gca; ax.GridColor = [0.2 0.2 0.2]; ax.GridAlpha = 0.1;
box off;
hold off;
end

function plot_perm_pvalues(R, pl)
% -log10(permutation p) per area
figure('Units', 'centimeters', 'Position', [1 8 20.4047 6], 'Color', 'w');
valid = ~isnan(R.Perm_p);
if ~any(valid), return; end
pvals = R.Perm_p(valid);
qvals = R.Perm_q(valid);
effs  = R.Perm_logOdds_U_vs_L(valid);
names = R.Area(valid);

[pvals_sorted, sidx] = sort_for_plot(pvals, pl.plot_seq, 'ascend');
qvals_sorted = qvals(sidx);
effs_sorted  = effs(sidx);
names_sorted = names(sidx);
neg_log_p = -log10(pvals_sorted);
nplot = numel(neg_log_p);
xpos = 1:nplot;

hold on;
for k = 1:nplot
    if qvals_sorted(k) < pl.QPlot
        if effs_sorted(k) > 0, col = pl.ColorU; else, col = pl.ColorL; end
    else
        col = pl.ColorNonSig;
    end
    bar(xpos(k), neg_log_p(k), 'FaceColor', col, 'EdgeColor', col * 0.6, 'LineWidth', 0.5);
end
yline(-log10(0.05), '--', 'p = 0.05', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xlim([0 nplot+1]);
set(gca, 'XTick', xpos, 'XTickLabel', names_sorted, 'FontSize', 10, ...
    'TickDir', 'out', 'XTickLabelRotation', 90);
ylabel('-log_{10}(p-value)', 'FontSize', 11);
title('Permutation Test: P-values by Area', 'FontSize', 12);
grid on;
ax = gca; ax.GridColor = [0.2 0.2 0.2]; ax.GridAlpha = 0.1;
ax.YGrid = 'on'; ax.XGrid = 'off';
box off;
hold off;
end

function plot_perm_boot_ci(R, pl)
% Permutation effects with Bonferroni BCa bootstrap CIs; colour by permutation q
figure('Units', 'centimeters', 'Position', [16.4888 1.5 20.4047 6.8580], 'Color', 'w');
valid = ~isnan(R.Perm_logOdds_U_vs_L) & ~isnan(R.Boot_CIlo_BCa_Bonf);
if ~any(valid)
    text(0.1, 0.5, 'No valid Permutation results to plot', 'FontSize', 11);
    return;
end
effs  = R.Perm_logOdds_U_vs_L(valid);
ci_lo = R.Boot_CIlo_BCa_Bonf(valid);
ci_hi = R.Boot_CIhi_BCa_Bonf(valid);
names = R.Area(valid);
abund_pct = 100 * (R.meanFracU(valid) + R.meanFracL(valid)) ./ 2;
qvals = R.Perm_q(valid);

[eff_sorted, sidx] = sort_for_plot(effs, pl.plot_seq, 'descend');
ci_lo = ci_lo(sidx); ci_hi = ci_hi(sidx); names = names(sidx);
abund_pct = abund_pct(sidx); q_sorted = qvals(sidx);

[sVals, legVals, legSizes] = abundance_marker_sizes(abund_pct);
cols = significance_colors(eff_sorted, q_sorted, pl);

draw_forest_columns(eff_sorted, ci_lo, ci_hi, names, sVals, cols, cols, 1, ...
    legVals, legSizes, 0.6, sprintf('Log-odds ratio (%s vs %s)', pl.GroupNames{:}), ...
    'Permutation p-values + Bootstrap CIs');
end

function plot_glmm_volcano(R, pl)
% log2 fold change vs -log10(GLMM p); coloured and labelled if GLMM q < 0.05
figure('Units', 'centimeters', 'Position', [2.4342 6.8157 13.5678 10.2658], 'Color', 'w');
ax1 = axes('Position', [0.1 0.15 0.60 0.78]); hold(ax1, 'on');
x = R.log2FC;
y = -log10(R.GLMM_p);
abundance = (R.meanFracU + R.meanFracL) ./ 2;
s_all = max(10, 1000 * abundance);
sig = R.GLMM_q < 0.05;
sig_pos = sig & x > 0;
sig_neg = sig & x <= 0;

scatter(ax1, x(~sig), y(~sig), s_all(~sig), 'MarkerFaceColor', pl.ColorNonSig, ...
    'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerFaceAlpha', 0.6);
scatter(ax1, x(sig_pos), y(sig_pos), s_all(sig_pos), 'MarkerFaceColor', pl.ColorU, ...
    'MarkerEdgeColor', pl.ColorU - 0.1, 'MarkerFaceAlpha', 1);
scatter(ax1, x(sig_neg), y(sig_neg), s_all(sig_neg), 'MarkerFaceColor', pl.ColorL, ...
    'MarkerEdgeColor', pl.ColorL - 0.1, 'MarkerFaceAlpha', 1);
xline(ax1, 0, '--', 'LineWidth', 1, 'Color', [0 0 0]);
xlabel(ax1, sprintf('log_2(%% %s / %% %s)', pl.GroupNames{:}), 'FontSize', 11);
ylabel(ax1, '-log_{10}(p-value)', 'FontSize', 11);
box(ax1, 'off'); grid(ax1, 'off'); set(ax1, 'TickDir', 'out', 'FontSize', 10);

labelIdx = find(sig & ~isnan(y));
textOffsetX = 0.03 * max(1, range(x(isfinite(x))));
for k = 1:numel(labelIdx)
    idx = labelIdx(k);
    text(ax1, x(idx) + textOffsetX, y(idx), R.Area{idx}, 'FontSize', 10);
end

% Colour legend
axPos = get(ax1, 'Position');
legendX = axPos(1) + axPos(3) + 0.02;
legendW = 0.18;
lgAx = axes('Position', [legendX, axPos(2) + 0.18, legendW, 0.12], 'Visible', 'off'); hold(lgAx, 'on');
scatter(lgAx, 0.23, 0.7, 80, 'o', 'MarkerFaceColor', pl.ColorU, 'MarkerEdgeColor', 'none');
text(lgAx, 0.38, 0.7, 'q < 0.05 (FDR)', 'FontSize', 9, 'VerticalAlignment', 'middle');
scatter(lgAx, 0.23, 0.25, 80, 'o', 'MarkerFaceColor', pl.ColorNonSig, 'MarkerEdgeColor', 'none');
text(lgAx, 0.38, 0.25, 'q >= 0.05', 'FontSize', 9, 'VerticalAlignment', 'middle');
xlim(lgAx, [0 1]); ylim(lgAx, [0 1]); axis(lgAx, 'off');

% Bubble-size legend
bubbleAx = axes('Position', [legendX, axPos(2), legendW, 0.18], 'Visible', 'off'); hold(bubbleAx, 'on');
ab_non = abundance(~isnan(abundance) & abundance >= 0);
if isempty(ab_non)
    pcts = [1 0.5 0.1] * 100;
else
    pcts = [max(ab_non), mean(ab_non), min(ab_non)] * 100;
end
svals = max(10, 1000 * (pcts / 100));
yspot = [0.75, 0.45, 0.15];
for k = 1:3
    scatter(bubbleAx, 0.23, yspot(k), svals(k), 'MarkerFaceColor', pl.ColorNonSig, ...
        'MarkerEdgeColor', [0.3 0.3 0.3]);
    text(bubbleAx, 0.38, yspot(k), sprintf('%.1f%%', pcts(k)), 'FontSize', 10, ...
        'VerticalAlignment', 'middle');
end
xlim(bubbleAx, [0 1]); ylim(bubbleAx, [0 1]); axis(bubbleAx, 'off');
hold(ax1, 'off');
end

function plot_glmm_forest_vertical(R, pl)
% GLMM effects with Bonferroni CIs, one row per area, sorted by effect
figure('Units', 'centimeters', 'Position', [12.2767 2.6882 11.7475 17.4202], 'Color', 'w');
valid = ~isnan(R.GLMM_logOdds_U_vs_L);
if ~any(valid)
    text(0.1, 0.5, 'No valid GLMM results to plot', 'FontSize', 11);
    return;
end
effs  = R.GLMM_logOdds_U_vs_L(valid);
ci_lo = R.GLMM_CIlo(valid);
ci_hi = R.GLMM_CIhi(valid);
names = R.Area(valid);
abund_pct = 100 * (R.meanFracU(valid) + R.meanFracL(valid)) ./ 2;

[eff_sorted, sidx] = sort(effs, 'descend');
ci_lo = ci_lo(sidx); ci_hi = ci_hi(sidx); names = names(sidx); abund_pct = abund_pct(sidx);
[sVals, legVals, legSizes] = abundance_marker_sizes(abund_pct);

neutralDot = [0.7 0.7 0.7];
nplot = numel(eff_sorted);
ypos = 1:nplot;
hold on;
for k = 1:nplot
    e = eff_sorted(k); lo = ci_lo(k); hi = ci_hi(k);
    if e > 0 && lo > 0
        col = pl.ColorU;
    elseif e <= 0 && hi < 0
        col = pl.ColorL;
    else
        col = neutralDot;
    end
    plot([lo hi], [ypos(k) ypos(k)], '-', 'LineWidth', 1.8, 'Color', col);
    scatter(e, ypos(k), sVals(k), 'o', 'MarkerFaceColor', col, ...
        'MarkerEdgeColor', col * 0.5, 'LineWidth', 0.5, 'MarkerFaceAlpha', 1);
end
plot([0 0], [0 nplot+1], '--k', 'LineWidth', 1);

pad = 0.5;
xlim([min(ci_lo) - pad, max(ci_hi) + pad]);
ylim([0.5 nplot + 0.5]);
set(gca, 'YTick', ypos, 'YTickLabel', names, 'FontSize', 10, 'TickDir', 'out');
ax = gca; ax.YDir = 'reverse';
box off;
set(gca, 'XGrid', 'on', 'YGrid', 'off');
ax.GridColor = [0.2 0.2 0.2];
ax.GridAlpha = 0.05;
xlabel(sprintf('log-odds (%s vs %s)', pl.GroupNames{:}), 'FontSize', 11);

% Marker-size legend
yspot = [23, 24, 25];
xlims = xlim;
for i = 1:3
    scatter(gca, xlims(2) / 2, yspot(i), legSizes(i), 'o', ...
        'MarkerFaceColor', neutralDot, 'MarkerEdgeColor', neutralDot);
    text(gca, xlims(2) / 2 + 0.5, yspot(i), sprintf('%.1f%%', legVals(i)), ...
        'FontSize', 11, 'VerticalAlignment', 'middle');
end
hold off;
end

function plot_glmm_forest(R, pl)
% GLMM effects with Bonferroni CIs, one column per area, coloured by effect size
figure('Units', 'centimeters', 'Position', [16.4888 8.4878 20.4047 6.8580], 'Color', 'w');
valid = ~isnan(R.GLMM_logOdds_U_vs_L);
if ~any(valid)
    text(0.1, 0.5, 'No valid GLMM results to plot', 'FontSize', 11);
    return;
end
effs  = R.GLMM_logOdds_U_vs_L(valid);
ci_lo = R.GLMM_CIlo(valid);
ci_hi = R.GLMM_CIhi(valid);
names = R.Area(valid);
abund_pct = 100 * (R.meanFracU(valid) + R.meanFracL(valid)) ./ 2;

[eff_sorted, sidx] = sort_for_plot(effs, pl.plot_seq, 'descend');
ci_lo = ci_lo(sidx); ci_hi = ci_hi(sidx); names = names(sidx); abund_pct = abund_pct(sidx);
[sVals, legVals, legSizes] = abundance_marker_sizes(abund_pct);

% Colour scaled by effect relative to the largest positive (negative) effect
maxPos = max([eps; eff_sorted(eff_sorted > 0)]);
maxNeg = max([eps; abs(eff_sorted(eff_sorted < 0))]);
nplot = numel(eff_sorted);
cols = repmat([0.7 0.7 0.7], nplot, 1);
for k = 1:nplot
    e = eff_sorted(k);
    if e > 0
        cols(k, :) = palette_color(min(1, e / maxPos), pl.PosPalette);
    elseif e <= 0
        cols(k, :) = palette_color(min(1, abs(e) / maxNeg), pl.NegPalette);
    end
end

draw_forest_columns(eff_sorted, ci_lo, ci_hi, names, sVals, cols, cols, 1, ...
    legVals, legSizes, 0.8, sprintf('log-odds (%s vs %s)', pl.GroupNames{:}), '');
end

function plot_glmm_boxplots(R, pl)
% Per-animal input fractions (U vs L) for areas with GLMM p < 0.05
sel = find(R.GLMM_p < 0.05);
if isempty(sel)
    qm = R.GLMM_q; qm(isnan(qm)) = Inf;
    [~, ord] = sort(qm, 'ascend');
    sel = ord(1:min(6, numel(ord)));
end
if isempty(sel), return; end

nU_anim = size(pl.fracU, 1);
nL_anim = size(pl.fracL, 1);
figure('Units', 'centimeters', 'Position', [12.7847 1.0795 15.8115 19.6003], 'Color', 'w');
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
boxW = 0.6;
for ii = 1:numel(sel)
    idx = sel(ii);
    ax = nexttile(t); hold(ax, 'on');
    dataU = pl.fracU(:, idx);
    dataL = pl.fracL(:, idx);

    rng(100 + ii);                                 % jitter of the points
    draw_box(ax, 1, dataU, pl.ColorU, boxW, nU_anim);
    draw_box(ax, 2, dataL, pl.ColorL, boxW, nL_anim);

    ylabel(ax, 'Input fraction', 'FontSize', 9);
    title(ax, sprintf('%s\n p=%.3g, q=%.3g', R.Area{idx}, R.GLMM_p(idx), R.GLMM_q(idx)), 'FontSize', 9);
    allVals = [dataU; dataL];
    rangeV = max(allVals) - min(allVals);
    if rangeV == 0, pad = 0.01; else, pad = 0.08 * rangeV; end
    ylim(ax, [min(allVals) - pad, max(allVals) + pad]);
    xlim(ax, [0.4 2.6]);
    set(ax, 'XTick', [1 2], 'XTickLabel', pl.GroupTicks, 'FontSize', 9, 'TickDir', 'out');
    box(ax, 'off'); grid(ax, 'off');
    hold(ax, 'off');
end
end

function draw_box(ax, xc, data, col, boxW, nAnim)
% Box (IQR), points, whiskers (1.5 IQR), median line and mean dot
q1 = prctile(data, 25); q3 = prctile(data, 75);
iqrV = q3 - q1;
whisk_lo = max(min(data), q1 - 1.5 * iqrV);
whisk_hi = min(max(data), q3 + 1.5 * iqrV);
patch([xc-boxW/2, xc+boxW/2, xc+boxW/2, xc-boxW/2], [q1, q1, q3, q3], col, ...
    'FaceAlpha', 0.35, 'EdgeColor', 'none', 'Parent', ax);
jitter = (rand(nAnim, 1) - 0.5) * 0.12;
scatter(ax, xc + jitter, data, 36, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.9);
plot(ax, [xc xc], [q3 whisk_hi], '-', 'Color', col * 0.8, 'LineWidth', 1.4);
plot(ax, [xc xc], [q1 whisk_lo], '-', 'Color', col * 0.8, 'LineWidth', 1.4);
plot(ax, [xc-boxW/2, xc+boxW/2], [median(data) median(data)], '-k', 'LineWidth', 2);
plot(ax, xc, mean(data), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
end

function plot_perm_inverted_ci(R, pl)
% Permutation effects with FCR-adjusted inverted permutation CIs
figure('Units', 'centimeters', 'Position', [16.4888 1.5 20.4047 6.8580], 'Color', 'w');
valid = ~isnan(R.Perm_logOdds_U_vs_L) & ~isnan(R.Perm_CIlo) & ~isnan(R.Perm_CIhi);
if ~any(valid)
    text(0.1, 0.5, 'No valid Permutation results to plot', 'FontSize', 11);
    return;
end
effs  = R.Perm_logOdds_U_vs_L(valid);
ci_lo = R.Perm_CIlo_FCR(valid);
ci_hi = R.Perm_CIhi_FCR(valid);
names = R.Area(valid);
abund_pct = 100 * (R.meanFracU(valid) + R.meanFracL(valid)) ./ 2;
qvals = R.Perm_q(valid);

[eff_sorted, sidx] = sort_for_plot(effs, pl.plot_seq, 'descend');
ci_lo = ci_lo(sidx); ci_hi = ci_hi(sidx); names = names(sidx);
abund_pct = abund_pct(sidx); q_sorted = qvals(sidx);

[sVals, legVals, legSizes] = abundance_marker_sizes(abund_pct);
cols = significance_colors(eff_sorted, q_sorted, pl);

draw_forest_columns(eff_sorted, ci_lo, ci_hi, names, sVals, cols, cols, 1, ...
    legVals, legSizes, 0.6, sprintf('Log-odds ratio (%s vs %s)', pl.GroupNames{:}), ...
    'Permutation test with inverted CI');
end
