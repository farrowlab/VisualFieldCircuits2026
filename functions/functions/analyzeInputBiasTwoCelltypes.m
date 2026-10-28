function Results = analyzeInputBiasTwoCelltypes(nu1, nl1, nu2, nl2, areaNames, varargin)
% ANALYZEINPUTBIASTWOCELLTYPES
% Compare upper vs lower inputs for two cell types and visualize.
%
% Usage:
% Results = analyzeInputBiasTwoCelltypes(nu1, nl1, nu2, nl2, areaNames, ...
%       'TotalU1',N_U1,'TotalL1',N_L1,'TotalU2',N_U2,'TotalL2',N_L2, ...
%       'ColorIntPos', [r g b], 'ColorIntNeg', [r g b], 'Plot', true);

weighted = false;  % true = pooled (weighted by trials), false = mean (unweighted, recommended)

% -------------------- Parse inputs --------------------
p = inputParser;
addRequired(p,'nu1',@isnumeric);
addRequired(p,'nl1',@isnumeric);
addRequired(p,'nu2',@isnumeric);
addRequired(p,'nl2',@isnumeric);
addRequired(p,'areaNames',@(x) iscellstr(x)||isstring(x)||ischar(x));
addParameter(p,'TotalU1',[],@(x) isnumeric(x) && (isscalar(x) || isvector(x)) || isempty(x));
addParameter(p,'TotalL1',[],@(x) isnumeric(x) && (isscalar(x) || isvector(x)) || isempty(x));
addParameter(p,'TotalU2',[],@(x) isnumeric(x) && (isscalar(x) || isvector(x)) || isempty(x));
addParameter(p,'TotalL2',[],@(x) isnumeric(x) && (isscalar(x) || isvector(x)) || isempty(x));
addParameter(p,'Alpha',1e-10,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'UseWeights','inverseTrials',@(x) any(validatestring(x,{'inverseTrials','none'})));
addParameter(p,'PosPalette', [0.8510 0.5000 0; 0.8808 0.6000 0.2000; 0.9106 0.7000 0.4000; 0.9404 0.8000 0.6000; 0.9702 0.9000 0.8000], @(x) validatePalette(x));
addParameter(p,'NegPalette', [0.2100 0.3900 0.7300; 0.3680 0.5120 0.7840; 0.5260 0.6340 0.8380; 0.6840 0.7560 0.8920; 0.8420 0.8780 0.9460], @(x) validatePalette(x));
addParameter(p,'ColorU1',[0.91 0.68 0.49],@(x) validateColor(x));
addParameter(p,'ColorL1',[0.43 0.61 0.77],@(x) validateColor(x));
addParameter(p,'ColorU2',[0.91 0.68 0.49],@(x) validateColor(x));
addParameter(p,'ColorL2',[0.43 0.61 0.77],@(x) validateColor(x));
addParameter(p,'ColorIntPos',[0.91 0.68 0.49],@(x) validateColor(x));
addParameter(p,'ColorIntNeg',[0.43 0.61 0.77],@(x) validateColor(x));
addParameter(p,'ColorNonSig',[0.78 0.78 0.78],@(x) validateColor(x));
addParameter(p,'Plot',true,@islogical);

parse(p,nu1,nl1,nu2,nl2,areaNames,varargin{:});
Alpha = p.Results.Alpha;
UseWeights = p.Results.UseWeights;
PosPalette = p.Results.PosPalette;
NegPalette = p.Results.NegPalette;
ColorU1 = p.Results.ColorU1;
ColorL1 = p.Results.ColorL1;
ColorU2 = p.Results.ColorU2;
ColorL2 = p.Results.ColorL2;
ColorIntPos = p.Results.ColorIntPos;
ColorIntNeg = p.Results.ColorIntNeg;
ColorNonSig = p.Results.ColorNonSig;
doPlot = p.Results.Plot;

% Print weighting mode
if weighted
    fprintf('Using WEIGHTED estimation (pooled proportions - biased toward high-trial animals)\n');
else
    fprintf('Using UNWEIGHTED estimation (mean of proportions - unbiased animal-level inference)\n');
end

% -------------------- Sizes & checks --------------------
areaNames = cellstr(areaNames(:));
nAreas = numel(areaNames);

[nU1, nAreasU1] = size(nu1);
[nL1, nAreasL1] = size(nl1);
[nU2, nAreasU2] = size(nu2);
[nL2, nAreasL2] = size(nl2);

if ~(nAreasU1==nAreasL1 && nAreasU1==nAreasU2 && nAreasU1==nAreasL2 && nAreasU1==nAreas)
    error('All count matrices must have same number of columns (areas). Provided areaNames length = %d, but matrices columns vary: %d %d %d %d', ...
        nAreas, nAreasU1, nAreasL1, nAreasU2, nAreasL2);
end

% -------------------- Totals per animal (Trials) --------------------
if isempty(p.Results.TotalU1)
    N_U1 = sum(nu1,2);
    warning('TotalU1 not provided: defaulting to sum(nu1,2).');
else
    if isscalar(p.Results.TotalU1), N_U1 = repmat(p.Results.TotalU1, nU1,1); else N_U1 = p.Results.TotalU1(:); end
    assert(numel(N_U1)==nU1,'TotalU1 length mismatch.');
end
if isempty(p.Results.TotalL1)
    N_L1 = sum(nl1,2);
    warning('TotalL1 not provided: defaulting to sum(nl1,2).');
else
    if isscalar(p.Results.TotalL1), N_L1 = repmat(p.Results.TotalL1, nL1,1); else N_L1 = p.Results.TotalL1(:); end
    assert(numel(N_L1)==nL1,'TotalL1 length mismatch.');
end

if isempty(p.Results.TotalU2)
    N_U2 = sum(nu2,2);
    warning('TotalU2 not provided: defaulting to sum(nu2,2).');
else
    if isscalar(p.Results.TotalU2), N_U2 = repmat(p.Results.TotalU2, nU2,1); else N_U2 = p.Results.TotalU2(:); end
    assert(numel(N_U2)==nU2,'TotalU2 length mismatch.');
end
if isempty(p.Results.TotalL2)
    N_L2 = sum(nl2,2);
    warning('TotalL2 not provided: defaulting to sum(nl2,2).');
else
    if isscalar(p.Results.TotalL2), N_L2 = repmat(p.Results.TotalL2, nL2,1); else N_L2 = p.Results.TotalL2(:); end
    assert(numel(N_L2)==nL2,'TotalL2 length mismatch.');
end

if any(N_U1==0) || any(N_L1==0) || any(N_U2==0) || any(N_L2==0)
    error('TotalU/TotalL must be > 0 for all animals.');
end

% -------------------- Fractions & log2 fold-change --------------------
fracU1 = bsxfun(@rdivide, nu1, N_U1);   % nU1 x nAreas
fracL1 = bsxfun(@rdivide, nl1, N_L1);
meanFracU_CT1 = mean(fracU1,1)';    % column
meanFracL_CT1 = mean(fracL1,1)';    % column
log2FC_CT1 = log2( (meanFracU_CT1 + Alpha) ./ (meanFracL_CT1 + Alpha) );

fracU2 = bsxfun(@rdivide, nu2, N_U2);
fracL2 = bsxfun(@rdivide, nl2, N_L2);
meanFracU_CT2 = mean(fracU2,1)';    % column
meanFracL_CT2 = mean(fracL2,1)';    % column
log2FC_CT2 = log2( (meanFracU_CT2 + Alpha) ./ (meanFracL_CT2 + Alpha) );
% 
% % -------------------- Initialize storage --------------------
% GLMM_eff_CT1 = nan(nAreas,1); GLMM_se_CT1 = nan(nAreas,1); GLMM_p_CT1 = nan(nAreas,1);
% GLMM_CIlo_CT1 = nan(nAreas,1); GLMM_CIhi_CT1 = nan(nAreas,1);
% ranksum_p_CT1 = nan(nAreas,1);
% 
% GLMM_eff_CT2 = nan(nAreas,1); GLMM_se_CT2 = nan(nAreas,1); GLMM_p_CT2 = nan(nAreas,1);
% GLMM_CIlo_CT2 = nan(nAreas,1); GLMM_CIhi_CT2 = nan(nAreas,1);
% ranksum_p_CT2 = nan(nAreas,1);
% 
% GLMM_int_eff = nan(nAreas,1); GLMM_int_se = nan(nAreas,1); GLMM_int_p = nan(nAreas,1);
% GLMM_int_CIlo = nan(nAreas,1); GLMM_int_CIhi = nan(nAreas,1);


% % -------------------- Per-area GLMM per celltype --------------------
% for i = 1:nAreas
%     % ============================================================
%     % CT1 GLMM (binomial on counts)
%     % ============================================================
%     try
%         succ_U1 = nu1(:,i);
%         succ_L1 = nl1(:,i);
% 
%         % Ensure trials are vectors matching rows
%         if isscalar(N_U1), trials_U1 = repmat(N_U1, nU1, 1); else trials_U1 = N_U1(:); end
%         if isscalar(N_L1), trials_L1 = repmat(N_L1, nL1, 1); else trials_L1 = N_L1(:); end
% 
%         successes = [succ_U1; succ_L1];
%         trials = [trials_U1; trials_L1];
%         grp = [repmat("U", nU1, 1); repmat("L", nL1, 1)];
%         Animal = categorical([ arrayfun(@(k) sprintf('CT1U_%d',k), (1:nU1)','UniformOutput',false); ...
%                                arrayfun(@(k) sprintf('CT1L_%d',k), (1:nL1)','UniformOutput',false) ]);
% 
%         tbl = table(successes, trials, categorical(grp,["U","L"]), Animal, ...
%                     'VariableNames', {'Successes','Trials','Group','Animal'});
% 
%         % Fit binomial GLME
%         glme = fitglme(tbl, 'Successes ~ Group', ...
%             'Distribution','Binomial', 'BinomialSize', tbl.Trials, 'FitMethod','Laplace','Verbose',0);
% 
%         % Check for overdispersion
%         p_hat = predict(glme);
%         mu_hat_counts = p_hat .* tbl.Trials;
%         residuals = successes - mu_hat_counts;
%         var_hat = tbl.Trials .* p_hat .* (1 - p_hat);
%         nObs = numel(successes);
%         nFix = height(glme.Coefficients);
%         pearson_disp = sum((residuals).^2 ./ max(var_hat, eps)) / max(nObs - nFix,1);
% 
%         if pearson_disp > 1.5
%             % Add observation-level random effect
%             Obs = categorical(string((1:nObs)'));
%             tbl.Obs = Obs;
%             glme = fitglme(tbl, 'Successes ~ Group + (1|Obs)', ...
%                 'Distribution','Binomial', 'BinomialSize', tbl.Trials, 'FitMethod','Laplace', 'Verbose',0);
%         end
% 
%         coefTbl = glme.Coefficients;
%         idx = strcmp(coefTbl.Name,'Group_L'); 
%         if ~any(idx), idx = contains(coefTbl.Name,'Group'); end
%         b = coefTbl.Estimate(idx); se = coefTbl.SE(idx); pval = coefTbl.pValue(idx);
%         GLMM_eff_CT1(i) = -b; GLMM_se_CT1(i) = se; GLMM_p_CT1(i) = pval;
%         GLMM_CIlo_CT1(i) = -b - 1.96*se; GLMM_CIhi_CT1(i) = -b + 1.96*se;
%     catch ME
%         warning('GLMM CT1 failed for area %d: %s', i, ME.message);
%     end
%     try ranksum_p_CT1(i) = ranksum(fracU1(:,i), fracL1(:,i)); catch, ranksum_p_CT1(i)=NaN; end
% 
%     % ============================================================
%     % CT2 GLMM (binomial on counts)
%     % ============================================================
%     try
%         succ_U2 = nu2(:,i);
%         succ_L2 = nl2(:,i);
% 
%         if isscalar(N_U2), trials_U2 = repmat(N_U2, nU2, 1); else trials_U2 = N_U2(:); end
%         if isscalar(N_L2), trials_L2 = repmat(N_L2, nL2, 1); else trials_L2 = N_L2(:); end
% 
%         successes = [succ_U2; succ_L2];
%         trials = [trials_U2; trials_L2];
%         grp = [repmat("U", nU2, 1); repmat("L", nL2, 1)];
%         Animal = categorical([ arrayfun(@(k) sprintf('CT2U_%d',k), (1:nU2)','UniformOutput',false); ...
%                                arrayfun(@(k) sprintf('CT2L_%d',k), (1:nL2)','UniformOutput',false) ]);
% 
%         tbl = table(successes, trials, categorical(grp,["U","L"]), Animal, ...
%                     'VariableNames', {'Successes','Trials','Group','Animal'});
% 
%         glme = fitglme(tbl, 'Successes ~ Group', ...
%             'Distribution','Binomial', 'BinomialSize', tbl.Trials, 'FitMethod','Laplace','Verbose',0);
% 
%         % Check for overdispersion
%         p_hat = predict(glme);
%         mu_hat_counts = p_hat .* tbl.Trials;
%         residuals = successes - mu_hat_counts;
%         var_hat = tbl.Trials .* p_hat .* (1 - p_hat);
%         nObs = numel(successes);
%         nFix = height(glme.Coefficients);
%         pearson_disp = sum((residuals).^2 ./ max(var_hat, eps)) / max(nObs - nFix,1);
% 
%         if pearson_disp > 1.5
%             Obs = categorical(string((1:nObs)'));
%             tbl.Obs = Obs;
%             glme = fitglme(tbl, 'Successes ~ Group + (1|Obs)', ...
%                 'Distribution','Binomial', 'BinomialSize', tbl.Trials, 'FitMethod','Laplace', 'Verbose',0);
%         end
% 
%         coefTbl = glme.Coefficients;
%         idx = strcmp(coefTbl.Name,'Group_L'); 
%         if ~any(idx), idx = contains(coefTbl.Name,'Group'); end
%         b = coefTbl.Estimate(idx); se = coefTbl.SE(idx); pval = coefTbl.pValue(idx);
%         GLMM_eff_CT2(i) = -b; GLMM_se_CT2(i) = se; GLMM_p_CT2(i) = pval;
%         GLMM_CIlo_CT2(i) = -b - 1.96*se; GLMM_CIhi_CT2(i) = -b + 1.96*se;
%     catch ME
%         warning('GLMM CT2 failed for area %d: %s', i, ME.message);
%     end
%     try ranksum_p_CT2(i) = ranksum(fracU2(:,i), fracL2(:,i)); catch, ranksum_p_CT2(i)=NaN; end
% end
% 
% % FDR per celltype
% valid1 = ~isnan(GLMM_p_CT1);
% q_GLMM_CT1 = nan(nAreas,1);
% if any(valid1), q_GLMM_CT1(valid1) = mafdr(GLMM_p_CT1(valid1),'BHFDR',true); end
% 
% valid2 = ~isnan(GLMM_p_CT2);
% q_GLMM_CT2 = nan(nAreas,1);
% if any(valid2), q_GLMM_CT2(valid2) = mafdr(GLMM_p_CT2(valid2),'BHFDR',true); end
% 
% % -------------------- Combined GLMM for interaction using counts --------------------
% for i = 1:nAreas
%     try
%         successes = [ nu1(:,i); nl1(:,i); nu2(:,i); nl2(:,i) ];
%         trials =    [ N_U1;    N_L1;    N_U2;    N_L2 ];
%         Group = categorical([ repmat("U",nU1,1); repmat("L",nL1,1); repmat("U",nU2,1); repmat("L",nL2,1) ], ["U","L"]);
%         CellType = categorical([ repmat("CT1",nU1+nL1,1); repmat("CT2",nU2+nL2,1) ], ["CT1","CT2"]);
%         Animal = categorical([ arrayfun(@(k) sprintf('CT1U_%d',k), (1:nU1)','UniformOutput',false); ...
%                                arrayfun(@(k) sprintf('CT1L_%d',k), (1:nL1)','UniformOutput',false); ...
%                                arrayfun(@(k) sprintf('CT2U_%d',k), (1:nU2)','UniformOutput',false); ...
%                                arrayfun(@(k) sprintf('CT2L_%d',k), (1:nL2)','UniformOutput',false) ]);
%         tbl = table(successes, trials, Group, CellType, Animal, 'VariableNames', {'Successes','Trials','Group','CellType','Animal'});
% 
%         if strcmp(UseWeights,'inverseTrials')
%             w = round(median(tbl.Trials)./tbl.Trials);
%             w(w<1) = 1;
%             glme = fitglme(tbl, 'Successes ~ Group*CellType', ...
%                 'Distribution','Binomial','BinomialSize',tbl.Trials,'Weights',w,'FitMethod','Laplace','Verbose',0);
%         else
%             glme = fitglme(tbl, 'Successes ~ Group*CellType', ...
%                 'Distribution','Binomial','BinomialSize',tbl.Trials,'FitMethod','Laplace','Verbose',0);
%         end
% 
%         % Check for overdispersion
%         p_hat = predict(glme);
%         mu_hat_counts = p_hat .* tbl.Trials;
%         residuals = successes - mu_hat_counts;
%         var_hat = tbl.Trials .* p_hat .* (1 - p_hat);
%         nObs = numel(successes);
%         nFix = height(glme.Coefficients);
%         pearson_disp = sum((residuals).^2 ./ max(var_hat, eps)) / max(nObs - nFix,1);
% 
%         if pearson_disp > 1.5
%             Obs = categorical(string((1:nObs)'));
%             tbl.Obs = Obs;
%             glme = fitglme(tbl, 'Successes ~ Group*CellType + (1|Obs)', ...
%                 'Distribution','Binomial', 'BinomialSize', tbl.Trials, 'FitMethod','Laplace', 'Verbose',0);
%         end
% 
% 
%         coefTbl = glme.Coefficients;
%         coefNames = coefTbl.Name;
% 
%         % Robustly find the interaction coefficient index
%         idx_inter = find( contains(coefNames,'Group') & contains(coefNames,'CellType'), 1 );
%         if isempty(idx_inter)
%             idx_inter = find( ~cellfun(@isempty, regexp(coefNames, '(?i)(group).*:.*(celltype)|(celltype).*:.*(group)')), 1 );
%         end
%         if isempty(idx_inter)
%             idx_inter = find(contains(coefNames,':') & (contains(coefNames,'Group') | contains(coefNames,'CellType')),1);
%         end
% 
%         if ~isempty(idx_inter)
%             b_inter = coefTbl.Estimate(idx_inter);
%             se_inter = coefTbl.SE(idx_inter);
%             p_inter = coefTbl.pValue(idx_inter);
%             GLMM_int_eff(i) = -b_inter;
%             GLMM_int_se(i) = se_inter;
%             GLMM_int_p(i) = p_inter;
%             GLMM_int_CIlo(i) = -b_inter - 1.96*se_inter;
%             GLMM_int_CIhi(i) = -b_inter + 1.96*se_inter;
%         else
%             warning('Interaction coefficient not found for area %d; leaving NaN.', i);
%         end
%     catch ME
%         warning('GLMM interaction failed for area %d: %s', i, ME.message);
%     end
% end
% 
% % FDR on interactions
% validI = ~isnan(GLMM_int_p);
% GLMM_int_q = nan(nAreas,1);
% if any(validI), GLMM_int_q(validI) = mafdr(GLMM_int_p(validI),'BHFDR',true); end
% 
% % -------------------- BONFERRONI CI ADJUSTMENTS --------------------
% alpha_fam = 0.05;
% alpha_bonf = alpha_fam / nAreas;
% z_bonf = -norminv(alpha_bonf / 2);
% 
% % CT1 Bonferroni
% GLMM_CIlo_CT1_bonf = nan(nAreas,1);
% GLMM_CIhi_CT1_bonf = nan(nAreas,1);
% GLMM_CIlo_CT1_bonf(valid1) = GLMM_eff_CT1(valid1) - z_bonf * GLMM_se_CT1(valid1);
% GLMM_CIhi_CT1_bonf(valid1) = GLMM_eff_CT1(valid1) + z_bonf * GLMM_se_CT1(valid1);
% 
% % CT2 Bonferroni
% GLMM_CIlo_CT2_bonf = nan(nAreas,1);
% GLMM_CIhi_CT2_bonf = nan(nAreas,1);
% GLMM_CIlo_CT2_bonf(valid2) = GLMM_eff_CT2(valid2) - z_bonf * GLMM_se_CT2(valid2);
% GLMM_CIhi_CT2_bonf(valid2) = GLMM_eff_CT2(valid2) + z_bonf * GLMM_se_CT2(valid2);
% 
% % Interaction Bonferroni
% GLMM_int_CIlo_bonf = nan(nAreas,1);
% GLMM_int_CIhi_bonf = nan(nAreas,1);
% GLMM_int_CIlo_bonf(validI) = GLMM_int_eff(validI) - z_bonf * GLMM_int_se(validI);
% GLMM_int_CIhi_bonf(validI) = GLMM_int_eff(validI) + z_bonf * GLMM_int_se(validI);
% 
% % -------------------- FCR CI ADJUSTMENTS --------------------
% % CT1 FCR
% m1 = sum(valid1);
% R1 = sum(q_GLMM_CT1(valid1) < 0.05);
% if R1 == 0
%     alpha_adj_CT1 = 0.05 / m1;
%     fprintf('CT1 FCR: No discoveries. Using Bonferroni.\n');
% else
%     alpha_adj_CT1 = 0.05 * (R1 / m1);
%     fprintf('CT1 FCR: %d discoveries. Adjusted alpha = %.5f.\n', R1, alpha_adj_CT1);
% end
% z_fcr_CT1 = -norminv(alpha_adj_CT1 / 2);
% 
% GLMM_CIlo_CT1_fcr = nan(nAreas,1);
% GLMM_CIhi_CT1_fcr = nan(nAreas,1);
% GLMM_CIlo_CT1_fcr(valid1) = GLMM_eff_CT1(valid1) - z_fcr_CT1 * GLMM_se_CT1(valid1);
% GLMM_CIhi_CT1_fcr(valid1) = GLMM_eff_CT1(valid1) + z_fcr_CT1 * GLMM_se_CT1(valid1);
% 
% % CT2 FCR
% m2 = sum(valid2);
% R2 = sum(q_GLMM_CT2(valid2) < 0.05);
% if R2 == 0
%     alpha_adj_CT2 = 0.05 / m2;
%     fprintf('CT2 FCR: No discoveries. Using Bonferroni.\n');
% else
%     alpha_adj_CT2 = 0.05 * (R2 / m2);
%     fprintf('CT2 FCR: %d discoveries. Adjusted alpha = %.5f.\n', R2, alpha_adj_CT2);
% end
% z_fcr_CT2 = -norminv(alpha_adj_CT2 / 2);
% 
% GLMM_CIlo_CT2_fcr = nan(nAreas,1);
% GLMM_CIhi_CT2_fcr = nan(nAreas,1);
% GLMM_CIlo_CT2_fcr(valid2) = GLMM_eff_CT2(valid2) - z_fcr_CT2 * GLMM_se_CT2(valid2);
% GLMM_CIhi_CT2_fcr(valid2) = GLMM_eff_CT2(valid2) + z_fcr_CT2 * GLMM_se_CT2(valid2);
% 
% % Interaction FCR
% mI = sum(validI);
% RI = sum(GLMM_int_q(validI) < 0.05);
% if RI == 0
%     alpha_adj_int = 0.05 / mI;
%     fprintf('Interaction FCR: No discoveries. Using Bonferroni.\n');
% else
%     alpha_adj_int = 0.05 * (RI / mI);
%     fprintf('Interaction FCR: %d discoveries. Adjusted alpha = %.5f.\n', RI, alpha_adj_int);
% end
% z_fcr_int = -norminv(alpha_adj_int / 2);
% 
% GLMM_int_CIlo_fcr = nan(nAreas,1);
% GLMM_int_CIhi_fcr = nan(nAreas,1);
% GLMM_int_CIlo_fcr(validI) = GLMM_int_eff(validI) - z_fcr_int * GLMM_int_se(validI);
% GLMM_int_CIhi_fcr(validI) = GLMM_int_eff(validI) + z_fcr_int * GLMM_int_se(validI);


% -------------------- BCa Bootstrap + Permutation for CT1, CT2, and Interaction --------------------
rng(42);
nBootstrap = 10000;
nPermutations = 10000;
eps_prop = 1e-10;

% ==================== CT1 Storage ====================
Boot_eff_CT1 = nan(nAreas, 1);
Boot_CIlo_bca_CT1 = nan(nAreas, 1);
Boot_CIhi_bca_CT1 = nan(nAreas, 1);
Boot_CIlo_CT1 = nan(nAreas, 1);
Boot_CIhi_CT1 = nan(nAreas, 1);
Boot_z0_CT1 = nan(nAreas, 1);
Boot_a_CT1 = nan(nAreas, 1);
Boot_distributions_CT1 = cell(nAreas, 1);

Perm_p_CT1 = nan(nAreas, 1);
Perm_eff_CT1 = nan(nAreas, 1);
Perm_distributions_CT1 = cell(nAreas, 1);

% ==================== CT2 Storage ====================
Boot_eff_CT2 = nan(nAreas, 1);
Boot_CIlo_bca_CT2 = nan(nAreas, 1);
Boot_CIhi_bca_CT2 = nan(nAreas, 1);
Boot_CIlo_CT2 = nan(nAreas, 1);
Boot_CIhi_CT2 = nan(nAreas, 1);
Boot_z0_CT2 = nan(nAreas, 1);
Boot_a_CT2 = nan(nAreas, 1);
Boot_distributions_CT2 = cell(nAreas, 1);

Perm_p_CT2 = nan(nAreas, 1);
Perm_eff_CT2 = nan(nAreas, 1);
Perm_distributions_CT2 = cell(nAreas, 1);

% ==================== Interaction Storage ====================
Boot_int_eff = nan(nAreas, 1);
Boot_int_CIlo_bca = nan(nAreas, 1);
Boot_int_CIhi_bca = nan(nAreas, 1);
Boot_int_CIlo = nan(nAreas, 1);
Boot_int_CIhi = nan(nAreas, 1);
Boot_int_z0 = nan(nAreas, 1);
Boot_int_a = nan(nAreas, 1);
Boot_int_distributions = cell(nAreas, 1);

Perm_int_p = nan(nAreas, 1);
Perm_int_eff = nan(nAreas, 1);
Perm_int_distributions = cell(nAreas, 1);

% ==================== Main Loop ====================
% ==================== Main Loop ====================
for i = 1:nAreas
    
    % Extract data for this area
    succ_U1 = nu1(:, i);
    succ_L1 = nl1(:, i);
    succ_U2 = nu2(:, i);
    succ_L2 = nl2(:, i);
    
    trials_U1 = N_U1(:);
    trials_L1 = N_L1(:);
    trials_U2 = N_U2(:);
    trials_L2 = N_L2(:);
    
    % Flags for valid data
    valid_CT1 = ~any(isnan(succ_U1)) && ~any(isnan(succ_L1)) && ~(sum(succ_U1) == 0 && sum(succ_L1) == 0);
    valid_CT2 = ~any(isnan(succ_U2)) && ~any(isnan(succ_L2)) && ~(sum(succ_U2) == 0 && sum(succ_L2) == 0);
    
    % ============================================================
    % CT1: BOOTSTRAP + PERMUTATION
    % ============================================================
    
    if valid_CT1
        try
            % --- Bootstrap ---
            bootstrap_CT1 = zeros(nBootstrap, 1);
            for b = 1:nBootstrap
                idx_U1 = randsample(nU1, nU1, true);
                idx_L1 = randsample(nL1, nL1, true);
                
                succ_U1_boot = succ_U1(idx_U1);
                succ_L1_boot = succ_L1(idx_L1);
                trials_U1_boot = trials_U1(idx_U1);
                trials_L1_boot = trials_L1(idx_L1);
                
                if weighted
                    % WEIGHTED: Pooled proportions
                    prop_U1_boot = sum(succ_U1_boot) / sum(trials_U1_boot);
                    prop_L1_boot = sum(succ_L1_boot) / sum(trials_L1_boot);
                else
                    % UNWEIGHTED: Mean of proportions
                    prop_U1_boot = mean(succ_U1_boot ./ trials_U1_boot);
                    prop_L1_boot = mean(succ_L1_boot ./ trials_L1_boot);
                end
                
                prop_U1_boot = max(eps_prop, min(1 - eps_prop, prop_U1_boot));
                prop_L1_boot = max(eps_prop, min(1 - eps_prop, prop_L1_boot));
                
                odds_U1_boot = prop_U1_boot / (1 - prop_U1_boot);
                odds_L1_boot = prop_L1_boot / (1 - prop_L1_boot);
                bootstrap_CT1(b) = log(odds_U1_boot / odds_L1_boot);
            end
            
            % Calculate original estimate
            if weighted
                % WEIGHTED
                prop_U1_orig = sum(succ_U1) / sum(trials_U1);
                prop_L1_orig = sum(succ_L1) / sum(trials_L1);
            else
                % UNWEIGHTED
                prop_U1_orig = mean(succ_U1 ./ trials_U1);
                prop_L1_orig = mean(succ_L1 ./ trials_L1);
            end
            
            odds_U1_orig = prop_U1_orig / (1 - prop_U1_orig);
            odds_L1_orig = prop_L1_orig / (1 - prop_L1_orig);
            theta_CT1 = log(odds_U1_orig / odds_L1_orig);
            
            Boot_eff_CT1(i) = theta_CT1;
            Boot_distributions_CT1{i} = bootstrap_CT1;
            
            % --- BCa: Bias correction ---
            prop_less = sum(bootstrap_CT1 < theta_CT1) / nBootstrap;
            if prop_less == 0
                z0 = -3;
            elseif prop_less == 1
                z0 = 3;
            else
                z0 = norminv(prop_less);
            end
            Boot_z0_CT1(i) = z0;
            
            % --- BCa: Acceleration (jackknife) ---
            jackknife_CT1 = zeros(nU1 + nL1, 1);
            
            for j = 1:nU1
                idx_jack = setdiff(1:nU1, j);
                
                if weighted
                    % WEIGHTED
                    prop_U1_jack = sum(succ_U1(idx_jack)) / sum(trials_U1(idx_jack));
                    prop_L1_jack = sum(succ_L1) / sum(trials_L1);
                else
                    % UNWEIGHTED
                    prop_U1_jack = mean(succ_U1(idx_jack) ./ trials_U1(idx_jack));
                    prop_L1_jack = mean(succ_L1 ./ trials_L1);
                end
                
                odds_U1_jack = prop_U1_jack / (1 - prop_U1_jack);
                odds_L1_jack = prop_L1_jack / (1 - prop_L1_jack);
                jackknife_CT1(j) = log(odds_U1_jack / odds_L1_jack);
            end
            
            for j = 1:nL1
                idx_jack = setdiff(1:nL1, j);
                
                if weighted
                    % WEIGHTED
                    prop_U1_jack = sum(succ_U1) / sum(trials_U1);
                    prop_L1_jack = sum(succ_L1(idx_jack)) / sum(trials_L1(idx_jack));
                else
                    % UNWEIGHTED
                    prop_U1_jack = mean(succ_U1 ./ trials_U1);
                    prop_L1_jack = mean(succ_L1(idx_jack) ./ trials_L1(idx_jack));
                end
                
                odds_U1_jack = prop_U1_jack / (1 - prop_U1_jack);
                odds_L1_jack = prop_L1_jack / (1 - prop_L1_jack);
                jackknife_CT1(nU1 + j) = log(odds_U1_jack / odds_L1_jack);
            end
            
            % Calculate acceleration
            jack_mean = mean(jackknife_CT1);
            numerator = sum((jack_mean - jackknife_CT1).^3);
            denominator = 6 * (sum((jack_mean - jackknife_CT1).^2))^(3/2);
            
            if denominator == 0 || isnan(denominator) || abs(numerator/denominator) > 1
                a = 0;
            else
                a = numerator / denominator;
            end
            Boot_a_CT1(i) = a;
            
            % --- BCa CI ---
            alpha_level = 0.05;
            z_alpha_lower = norminv(alpha_level / 2);
            z_alpha_upper = norminv(1 - alpha_level / 2);
            
            alpha_lower_adj = normcdf(z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)));
            alpha_upper_adj = normcdf(z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)));
            
            alpha_lower_adj = max(0.001, min(0.999, alpha_lower_adj));
            alpha_upper_adj = max(0.001, min(0.999, alpha_upper_adj));
            
            bca_CI = prctile(bootstrap_CT1, [alpha_lower_adj * 100, alpha_upper_adj * 100]);
            Boot_CIlo_bca_CT1(i) = bca_CI(1);
            Boot_CIhi_bca_CT1(i) = bca_CI(2);
            
            percentile_CI = prctile(bootstrap_CT1, [2.5, 97.5]);
            Boot_CIlo_CT1(i) = percentile_CI(1);
            Boot_CIhi_CT1(i) = percentile_CI(2);
            
            % --- Permutation test ---
            all_succ_CT1 = [succ_U1; succ_L1];
            all_trials_CT1 = [trials_U1; trials_L1];
            n_total_CT1 = nU1 + nL1;
            
            Perm_eff_CT1(i) = theta_CT1;
            
            perm_CT1 = zeros(nPermutations, 1);
            for perm = 1:nPermutations
                perm_idx = randperm(n_total_CT1);
                
                perm_succ_U1 = all_succ_CT1(perm_idx(1:nU1));
                perm_succ_L1 = all_succ_CT1(perm_idx(nU1+1:end));
                perm_trials_U1 = all_trials_CT1(perm_idx(1:nU1));
                perm_trials_L1 = all_trials_CT1(perm_idx(nU1+1:end));
                
                if weighted
                    % WEIGHTED
                    prop_U1_perm = sum(perm_succ_U1) / sum(perm_trials_U1);
                    prop_L1_perm = sum(perm_succ_L1) / sum(perm_trials_L1);
                else
                    % UNWEIGHTED
                    prop_U1_perm = mean(perm_succ_U1 ./ perm_trials_U1);
                    prop_L1_perm = mean(perm_succ_L1 ./ perm_trials_L1);
                end
                
                prop_U1_perm = max(eps_prop, min(1 - eps_prop, prop_U1_perm));
                prop_L1_perm = max(eps_prop, min(1 - eps_prop, prop_L1_perm));
                
                odds_U1_perm = prop_U1_perm / (1 - prop_U1_perm);
                odds_L1_perm = prop_L1_perm / (1 - prop_L1_perm);
                perm_CT1(perm) = log(odds_U1_perm / odds_L1_perm);
            end
            
            Perm_p_CT1(i) = sum(abs(perm_CT1) >= abs(theta_CT1)) / nPermutations;
            Perm_distributions_CT1{i} = perm_CT1;
            
        catch ME
            warning('CT1 bootstrap/permutation failed for area %d: %s', i, ME.message);
        end
    end
    
    % ============================================================
    % CT2: BOOTSTRAP + PERMUTATION
    % ============================================================
    
    if valid_CT2
        try
            % --- Bootstrap ---
            bootstrap_CT2 = zeros(nBootstrap, 1);
            for b = 1:nBootstrap
                idx_U2 = randsample(nU2, nU2, true);
                idx_L2 = randsample(nL2, nL2, true);
                
                succ_U2_boot = succ_U2(idx_U2);
                succ_L2_boot = succ_L2(idx_L2);
                trials_U2_boot = trials_U2(idx_U2);
                trials_L2_boot = trials_L2(idx_L2);
                
                if weighted
                    % WEIGHTED: Pooled proportions
                    prop_U2_boot = sum(succ_U2_boot) / sum(trials_U2_boot);
                    prop_L2_boot = sum(succ_L2_boot) / sum(trials_L2_boot);
                else
                    % UNWEIGHTED: Mean of proportions
                    prop_U2_boot = mean(succ_U2_boot ./ trials_U2_boot);
                    prop_L2_boot = mean(succ_L2_boot ./ trials_L2_boot);
                end
                
                prop_U2_boot = max(eps_prop, min(1 - eps_prop, prop_U2_boot));
                prop_L2_boot = max(eps_prop, min(1 - eps_prop, prop_L2_boot));
                
                odds_U2_boot = prop_U2_boot / (1 - prop_U2_boot);
                odds_L2_boot = prop_L2_boot / (1 - prop_L2_boot);
                bootstrap_CT2(b) = log(odds_U2_boot / odds_L2_boot);
            end
            
            % Calculate original estimate
            if weighted
                % WEIGHTED
                prop_U2_orig = sum(succ_U2) / sum(trials_U2);
                prop_L2_orig = sum(succ_L2) / sum(trials_L2);
            else
                % UNWEIGHTED
                prop_U2_orig = mean(succ_U2 ./ trials_U2);
                prop_L2_orig = mean(succ_L2 ./ trials_L2);
            end
            
            odds_U2_orig = prop_U2_orig / (1 - prop_U2_orig);
            odds_L2_orig = prop_L2_orig / (1 - prop_L2_orig);
            theta_CT2 = log(odds_U2_orig / odds_L2_orig);
            
            Boot_eff_CT2(i) = theta_CT2;
            Boot_distributions_CT2{i} = bootstrap_CT2;
            
            % --- BCa: Bias correction ---
            prop_less = sum(bootstrap_CT2 < theta_CT2) / nBootstrap;
            if prop_less == 0
                z0 = -3;
            elseif prop_less == 1
                z0 = 3;
            else
                z0 = norminv(prop_less);
            end
            Boot_z0_CT2(i) = z0;
            
            % --- BCa: Acceleration (jackknife) ---
            jackknife_CT2 = zeros(nU2 + nL2, 1);
            
            for j = 1:nU2
                idx_jack = setdiff(1:nU2, j);
                
                if weighted
                    % WEIGHTED
                    prop_U2_jack = sum(succ_U2(idx_jack)) / sum(trials_U2(idx_jack));
                    prop_L2_jack = sum(succ_L2) / sum(trials_L2);
                else
                    % UNWEIGHTED
                    prop_U2_jack = mean(succ_U2(idx_jack) ./ trials_U2(idx_jack));
                    prop_L2_jack = mean(succ_L2 ./ trials_L2);
                end
                
                odds_U2_jack = prop_U2_jack / (1 - prop_U2_jack);
                odds_L2_jack = prop_L2_jack / (1 - prop_L2_jack);
                jackknife_CT2(j) = log(odds_U2_jack / odds_L2_jack);
            end
            
            for j = 1:nL2
                idx_jack = setdiff(1:nL2, j);
                
                if weighted
                    % WEIGHTED
                    prop_U2_jack = sum(succ_U2) / sum(trials_U2);
                    prop_L2_jack = sum(succ_L2(idx_jack)) / sum(trials_L2(idx_jack));
                else
                    % UNWEIGHTED
                    prop_U2_jack = mean(succ_U2 ./ trials_U2);
                    prop_L2_jack = mean(succ_L2(idx_jack) ./ trials_L2(idx_jack));
                end
                
                odds_U2_jack = prop_U2_jack / (1 - prop_U2_jack);
                odds_L2_jack = prop_L2_jack / (1 - prop_L2_jack);
                jackknife_CT2(nU2 + j) = log(odds_U2_jack / odds_L2_jack);
            end
            
            % Calculate acceleration
            jack_mean = mean(jackknife_CT2);
            numerator = sum((jack_mean - jackknife_CT2).^3);
            denominator = 6 * (sum((jack_mean - jackknife_CT2).^2))^(3/2);
            
            if denominator == 0 || isnan(denominator) || abs(numerator/denominator) > 1
                a = 0;
            else
                a = numerator / denominator;
            end
            Boot_a_CT2(i) = a;
            
            % --- BCa CI ---
            alpha_level = 0.05;
            z_alpha_lower = norminv(alpha_level / 2);
            z_alpha_upper = norminv(1 - alpha_level / 2);
            
            alpha_lower_adj = normcdf(z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)));
            alpha_upper_adj = normcdf(z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)));
            
            alpha_lower_adj = max(0.001, min(0.999, alpha_lower_adj));
            alpha_upper_adj = max(0.001, min(0.999, alpha_upper_adj));
            
            bca_CI = prctile(bootstrap_CT2, [alpha_lower_adj * 100, alpha_upper_adj * 100]);
            Boot_CIlo_bca_CT2(i) = bca_CI(1);
            Boot_CIhi_bca_CT2(i) = bca_CI(2);
            
            percentile_CI = prctile(bootstrap_CT2, [2.5, 97.5]);
            Boot_CIlo_CT2(i) = percentile_CI(1);
            Boot_CIhi_CT2(i) = percentile_CI(2);
            
            % --- Permutation test ---
            all_succ_CT2 = [succ_U2; succ_L2];
            all_trials_CT2 = [trials_U2; trials_L2];
            n_total_CT2 = nU2 + nL2;
            
            Perm_eff_CT2(i) = theta_CT2;
            
            perm_CT2 = zeros(nPermutations, 1);
            for perm = 1:nPermutations
                perm_idx = randperm(n_total_CT2);
                
                perm_succ_U2 = all_succ_CT2(perm_idx(1:nU2));
                perm_succ_L2 = all_succ_CT2(perm_idx(nU2+1:end));
                perm_trials_U2 = all_trials_CT2(perm_idx(1:nU2));
                perm_trials_L2 = all_trials_CT2(perm_idx(nU2+1:end));
                
                if weighted
                    % WEIGHTED
                    prop_U2_perm = sum(perm_succ_U2) / sum(perm_trials_U2);
                    prop_L2_perm = sum(perm_succ_L2) / sum(perm_trials_L2);
                else
                    % UNWEIGHTED
                    prop_U2_perm = mean(perm_succ_U2 ./ perm_trials_U2);
                    prop_L2_perm = mean(perm_succ_L2 ./ perm_trials_L2);
                end
                
                prop_U2_perm = max(eps_prop, min(1 - eps_prop, prop_U2_perm));
                prop_L2_perm = max(eps_prop, min(1 - eps_prop, prop_L2_perm));
                
                odds_U2_perm = prop_U2_perm / (1 - prop_U2_perm);
                odds_L2_perm = prop_L2_perm / (1 - prop_L2_perm);
                perm_CT2(perm) = log(odds_U2_perm / odds_L2_perm);
            end
            
            Perm_p_CT2(i) = sum(abs(perm_CT2) >= abs(theta_CT2)) / nPermutations;
            Perm_distributions_CT2{i} = perm_CT2;
            
        catch ME
            warning('CT2 bootstrap/permutation failed for area %d: %s', i, ME.message);
        end
    end
    
    % ============================================================
    % INTERACTION: BOOTSTRAP + PERMUTATION
    % ============================================================
    
    if valid_CT1 && valid_CT2
        try
            % --- Bootstrap ---
            bootstrap_int = zeros(nBootstrap, 1);
            for b = 1:nBootstrap
                idx_U1 = randsample(nU1, nU1, true);
                idx_L1 = randsample(nL1, nL1, true);
                idx_U2 = randsample(nU2, nU2, true);
                idx_L2 = randsample(nL2, nL2, true);
                
                succ_U1_boot = succ_U1(idx_U1);
                succ_L1_boot = succ_L1(idx_L1);
                trials_U1_boot = trials_U1(idx_U1);
                trials_L1_boot = trials_L1(idx_L1);
                
                succ_U2_boot = succ_U2(idx_U2);
                succ_L2_boot = succ_L2(idx_L2);
                trials_U2_boot = trials_U2(idx_U2);
                trials_L2_boot = trials_L2(idx_L2);
                
                if weighted
                    % WEIGHTED
                    prop_U1_boot = sum(succ_U1_boot) / sum(trials_U1_boot);
                    prop_L1_boot = sum(succ_L1_boot) / sum(trials_L1_boot);
                    prop_U2_boot = sum(succ_U2_boot) / sum(trials_U2_boot);
                    prop_L2_boot = sum(succ_L2_boot) / sum(trials_L2_boot);
                else
                    % UNWEIGHTED
                    prop_U1_boot = mean(succ_U1_boot ./ trials_U1_boot);
                    prop_L1_boot = mean(succ_L1_boot ./ trials_L1_boot);
                    prop_U2_boot = mean(succ_U2_boot ./ trials_U2_boot);
                    prop_L2_boot = mean(succ_L2_boot ./ trials_L2_boot);
                end
                
                prop_U1_boot = max(eps_prop, min(1 - eps_prop, prop_U1_boot));
                prop_L1_boot = max(eps_prop, min(1 - eps_prop, prop_L1_boot));
                prop_U2_boot = max(eps_prop, min(1 - eps_prop, prop_U2_boot));
                prop_L2_boot = max(eps_prop, min(1 - eps_prop, prop_L2_boot));
                
                odds_U1_boot = prop_U1_boot / (1 - prop_U1_boot);
                odds_L1_boot = prop_L1_boot / (1 - prop_L1_boot);
                logOR_CT1_boot = log(odds_U1_boot / odds_L1_boot);
                
                odds_U2_boot = prop_U2_boot / (1 - prop_U2_boot);
                odds_L2_boot = prop_L2_boot / (1 - prop_L2_boot);
                logOR_CT2_boot = log(odds_U2_boot / odds_L2_boot);
                
                bootstrap_int(b) = logOR_CT1_boot - logOR_CT2_boot;
            end
            
            % Calculate original estimate
            theta_int = theta_CT1 - theta_CT2;
            
            Boot_int_eff(i) = theta_int;
            Boot_int_distributions{i} = bootstrap_int;
            
            % --- BCa: Bias correction ---
            prop_less = sum(bootstrap_int < theta_int) / nBootstrap;
            if prop_less == 0
                z0 = -3;
            elseif prop_less == 1
                z0 = 3;
            else
                z0 = norminv(prop_less);
            end
            Boot_int_z0(i) = z0;
            
            % --- BCa: Acceleration (jackknife) ---
            n_total_animals = nU1 + nL1 + nU2 + nL2;
            jackknife_int = zeros(n_total_animals, 1);
            jack_idx = 0;
            
            % Use original proportions
            if weighted
                % WEIGHTED
                prop_U1_orig = sum(succ_U1) / sum(trials_U1);
                prop_L1_orig = sum(succ_L1) / sum(trials_L1);
                prop_U2_orig = sum(succ_U2) / sum(trials_U2);
                prop_L2_orig = sum(succ_L2) / sum(trials_L2);
            else
                % UNWEIGHTED
                prop_U1_orig = mean(succ_U1 ./ trials_U1);
                prop_L1_orig = mean(succ_L1 ./ trials_L1);
                prop_U2_orig = mean(succ_U2 ./ trials_U2);
                prop_L2_orig = mean(succ_L2 ./ trials_L2);
            end
            
            odds_U1_orig = prop_U1_orig / (1 - prop_U1_orig);
            odds_L1_orig = prop_L1_orig / (1 - prop_L1_orig);
            odds_U2_orig = prop_U2_orig / (1 - prop_U2_orig);
            odds_L2_orig = prop_L2_orig / (1 - prop_L2_orig);
            
            logOR_CT1_orig = log(odds_U1_orig / odds_L1_orig);
            logOR_CT2_orig = log(odds_U2_orig / odds_L2_orig);
            
            % Jackknife U1
            for j = 1:nU1
                jack_idx = jack_idx + 1;
                idx_jack = setdiff(1:nU1, j);
                
                if weighted
                    prop_U1_jack = sum(succ_U1(idx_jack)) / sum(trials_U1(idx_jack));
                else
                    prop_U1_jack = mean(succ_U1(idx_jack) ./ trials_U1(idx_jack));
                end
                
                odds_U1_jack = prop_U1_jack / (1 - prop_U1_jack);
                logOR_CT1_jack = log(odds_U1_jack / odds_L1_orig);
                jackknife_int(jack_idx) = logOR_CT1_jack - logOR_CT2_orig;
            end
            
            % Jackknife L1
            for j = 1:nL1
                jack_idx = jack_idx + 1;
                idx_jack = setdiff(1:nL1, j);
                
                if weighted
                    prop_L1_jack = sum(succ_L1(idx_jack)) / sum(trials_L1(idx_jack));
                else
                    prop_L1_jack = mean(succ_L1(idx_jack) ./ trials_L1(idx_jack));
                end
                
                odds_L1_jack = prop_L1_jack / (1 - prop_L1_jack);
                logOR_CT1_jack = log(odds_U1_orig / odds_L1_jack);
                jackknife_int(jack_idx) = logOR_CT1_jack - logOR_CT2_orig;
            end
            
            % Jackknife U2
            for j = 1:nU2
                jack_idx = jack_idx + 1;
                idx_jack = setdiff(1:nU2, j);
                
                if weighted
                    prop_U2_jack = sum(succ_U2(idx_jack)) / sum(trials_U2(idx_jack));
                else
                    prop_U2_jack = mean(succ_U2(idx_jack) ./ trials_U2(idx_jack));
                end
                
                odds_U2_jack = prop_U2_jack / (1 - prop_U2_jack);
                logOR_CT2_jack = log(odds_U2_jack / odds_L2_orig);
                jackknife_int(jack_idx) = logOR_CT1_orig - logOR_CT2_jack;
            end
            
            % Jackknife L2
            for j = 1:nL2
                jack_idx = jack_idx + 1;
                idx_jack = setdiff(1:nL2, j);
                
                if weighted
                    prop_L2_jack = sum(succ_L2(idx_jack)) / sum(trials_L2(idx_jack));
                else
                    prop_L2_jack = mean(succ_L2(idx_jack) ./ trials_L2(idx_jack));
                end
                
                odds_L2_jack = prop_L2_jack / (1 - prop_L2_jack);
                logOR_CT2_jack = log(odds_U2_orig / odds_L2_jack);
                jackknife_int(jack_idx) = logOR_CT1_orig - logOR_CT2_jack;
            end
            
            % Calculate acceleration
            jack_mean = mean(jackknife_int);
            numerator = sum((jack_mean - jackknife_int).^3);
            denominator = 6 * (sum((jack_mean - jackknife_int).^2))^(3/2);
            
            if denominator == 0 || isnan(denominator) || abs(numerator/denominator) > 1
                a = 0;
            else
                a = numerator / denominator;
            end
            Boot_int_a(i) = a;
            
            % --- BCa CI ---
            alpha_level = 0.05;
            z_alpha_lower = norminv(alpha_level / 2);
            z_alpha_upper = norminv(1 - alpha_level / 2);
            
            alpha_lower_adj = normcdf(z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)));
            alpha_upper_adj = normcdf(z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)));
            
            alpha_lower_adj = max(0.001, min(0.999, alpha_lower_adj));
            alpha_upper_adj = max(0.001, min(0.999, alpha_upper_adj));
            
            bca_CI = prctile(bootstrap_int, [alpha_lower_adj * 100, alpha_upper_adj * 100]);
            Boot_int_CIlo_bca(i) = bca_CI(1);
            Boot_int_CIhi_bca(i) = bca_CI(2);
            
            percentile_CI = prctile(bootstrap_int, [2.5, 97.5]);
            Boot_int_CIlo(i) = percentile_CI(1);
            Boot_int_CIhi(i) = percentile_CI(2);
            
            % --- Permutation test for interaction ---
            all_succ_U = [succ_U1; succ_U2];
            all_trials_U = [trials_U1; trials_U2];
            all_succ_L = [succ_L1; succ_L2];
            all_trials_L = [trials_L1; trials_L2];
            
            Perm_int_eff(i) = theta_int;
            
            perm_int = zeros(nPermutations, 1);
            for perm = 1:nPermutations
                % Shuffle within Upper group
                perm_idx_U = randperm(nU1 + nU2);
                perm_succ_U1 = all_succ_U(perm_idx_U(1:nU1));
                perm_succ_U2 = all_succ_U(perm_idx_U(nU1+1:end));
                perm_trials_U1 = all_trials_U(perm_idx_U(1:nU1));
                perm_trials_U2 = all_trials_U(perm_idx_U(nU1+1:end));
                
                % Shuffle within Lower group
                perm_idx_L = randperm(nL1 + nL2);
                perm_succ_L1 = all_succ_L(perm_idx_L(1:nL1));
                perm_succ_L2 = all_succ_L(perm_idx_L(nL1+1:end));
                perm_trials_L1 = all_trials_L(perm_idx_L(1:nL1));
                perm_trials_L2 = all_trials_L(perm_idx_L(nL1+1:end));
                
                if weighted
                    % WEIGHTED
                    prop_U1_perm = sum(perm_succ_U1) / sum(perm_trials_U1);
                    prop_L1_perm = sum(perm_succ_L1) / sum(perm_trials_L1);
                    prop_U2_perm = sum(perm_succ_U2) / sum(perm_trials_U2);
                    prop_L2_perm = sum(perm_succ_L2) / sum(perm_trials_L2);
                else
                    % UNWEIGHTED
                    prop_U1_perm = mean(perm_succ_U1 ./ perm_trials_U1);
                    prop_L1_perm = mean(perm_succ_L1 ./ perm_trials_L1);
                    prop_U2_perm = mean(perm_succ_U2 ./ perm_trials_U2);
                    prop_L2_perm = mean(perm_succ_L2 ./ perm_trials_L2);
                end
                
                prop_U1_perm = max(eps_prop, min(1 - eps_prop, prop_U1_perm));
                prop_L1_perm = max(eps_prop, min(1 - eps_prop, prop_L1_perm));
                prop_U2_perm = max(eps_prop, min(1 - eps_prop, prop_U2_perm));
                prop_L2_perm = max(eps_prop, min(1 - eps_prop, prop_L2_perm));
                
                odds_U1_perm = prop_U1_perm / (1 - prop_U1_perm);
                odds_L1_perm = prop_L1_perm / (1 - prop_L1_perm);
                logOR_CT1_perm = log(odds_U1_perm / odds_L1_perm);
                
                odds_U2_perm = prop_U2_perm / (1 - prop_U2_perm);
                odds_L2_perm = prop_L2_perm / (1 - prop_L2_perm);
                logOR_CT2_perm = log(odds_U2_perm / odds_L2_perm);
                
                perm_int(perm) = logOR_CT1_perm - logOR_CT2_perm;
            end
            
            Perm_int_p(i) = sum(abs(perm_int) >= abs(theta_int)) / nPermutations;
            Perm_int_distributions{i} = perm_int;
            
        catch ME
            warning('Interaction bootstrap/permutation failed for area %d: %s', i, ME.message);
        end
    end
end

% ==================== FDR Correction for Permutation P-values ====================
% CT1
valid_CT1 = ~isnan(Perm_p_CT1);
Perm_q_CT1 = nan(nAreas, 1);
if any(valid_CT1)
    Perm_q_CT1(valid_CT1) = mafdr(Perm_p_CT1(valid_CT1), 'BHFDR', true);
end

% CT2
valid_CT2 = ~isnan(Perm_p_CT2);
Perm_q_CT2 = nan(nAreas, 1);
if any(valid_CT2)
    Perm_q_CT2(valid_CT2) = mafdr(Perm_p_CT2(valid_CT2), 'BHFDR', true);
end

% Interaction
valid_int = ~isnan(Perm_int_p);
Perm_int_q = nan(nAreas, 1);
if any(valid_int)
    Perm_int_q(valid_int) = mafdr(Perm_int_p(valid_int), 'BHFDR', true);
end

% ==================== FCR Adjustment for Bootstrap CIs ====================
% NOTE: FCR adjustment is based on PERMUTATION p-values, not bootstrap p-values

% --- CT1 FCR ---
validIdx_perm_CT1 = ~isnan(Perm_p_CT1);
m_CT1 = sum(validIdx_perm_CT1);

% Calculate q-values based on permutation p-values
qvals_CT1 = nan(nAreas, 1);
if m_CT1 > 0
    qvals_CT1(validIdx_perm_CT1) = mafdr(Perm_p_CT1(validIdx_perm_CT1), 'BHFDR', true);
end
Perm_q_CT1 = qvals_CT1;

% Count discoveries
R_CT1 = sum(qvals_CT1(validIdx_perm_CT1) < 0.045);

% Calculate FCR-adjusted alpha
if R_CT1 == 0
    alpha_adj_CT1 = 0.05 / m_CT1;
    fprintf('CT1 FCR: No discoveries. Using Bonferroni.\n');
else
    alpha_adj_CT1 = 0.05 * (R_CT1 / m_CT1);
    fprintf('CT1 FCR: %d discoveries. Adjusted alpha = %.5f.\n', R_CT1, alpha_adj_CT1);
end

% Adjust bootstrap BCa CIs using FCR
Boot_CIlo_bca_fcr_CT1 = nan(nAreas, 1);
Boot_CIhi_bca_fcr_CT1 = nan(nAreas, 1);
Boot_CIlo_bca_bonf_CT1 = nan(nAreas, 1);
Boot_CIhi_bca_bonf_CT1 = nan(nAreas, 1);
alpha_bonf = 0.05 / nAreas;

for i = 1:nAreas
    if ~isempty(Boot_distributions_CT1{i}) && ~isnan(Boot_z0_CT1(i))
        z0 = Boot_z0_CT1(i);
        a = Boot_a_CT1(i);
        
        % FCR adjustment
        z_alpha_lower_fcr = norminv(alpha_adj_CT1 / 2);
        z_alpha_upper_fcr = norminv(1 - alpha_adj_CT1 / 2);
        
        alpha_lower_adj_fcr = normcdf(z0 + (z0 + z_alpha_lower_fcr) / (1 - a * (z0 + z_alpha_lower_fcr)));
        alpha_upper_adj_fcr = normcdf(z0 + (z0 + z_alpha_upper_fcr) / (1 - a * (z0 + z_alpha_upper_fcr)));
        
        alpha_lower_adj_fcr = max(0.001, min(0.999, alpha_lower_adj_fcr));
        alpha_upper_adj_fcr = max(0.001, min(0.999, alpha_upper_adj_fcr));
        
        ci_fcr = prctile(Boot_distributions_CT1{i}, [alpha_lower_adj_fcr * 100, alpha_upper_adj_fcr * 100]);
        Boot_CIlo_bca_fcr_CT1(i) = ci_fcr(1);
        Boot_CIhi_bca_fcr_CT1(i) = ci_fcr(2);
        
        % Bonferroni adjustment
        z_alpha_lower_bonf = norminv(alpha_bonf / 2);
        z_alpha_upper_bonf = norminv(1 - alpha_bonf / 2);
        
        alpha_lower_adj_bonf = normcdf(z0 + (z0 + z_alpha_lower_bonf) / (1 - a * (z0 + z_alpha_lower_bonf)));
        alpha_upper_adj_bonf = normcdf(z0 + (z0 + z_alpha_upper_bonf) / (1 - a * (z0 + z_alpha_upper_bonf)));
        
        alpha_lower_adj_bonf = max(0.001, min(0.999, alpha_lower_adj_bonf));
        alpha_upper_adj_bonf = max(0.001, min(0.999, alpha_upper_adj_bonf));
        
        ci_bonf = prctile(Boot_distributions_CT1{i}, [alpha_lower_adj_bonf * 100, alpha_upper_adj_bonf * 100]);
        Boot_CIlo_bca_bonf_CT1(i) = ci_bonf(1);
        Boot_CIhi_bca_bonf_CT1(i) = ci_bonf(2);
    end
end

% --- CT2 FCR ---
validIdx_perm_CT2 = ~isnan(Perm_p_CT2);
m_CT2 = sum(validIdx_perm_CT2);

% Calculate q-values based on permutation p-values
qvals_CT2 = nan(nAreas, 1);
if m_CT2 > 0
    qvals_CT2(validIdx_perm_CT2) = mafdr(Perm_p_CT2(validIdx_perm_CT2), 'BHFDR', true);
end
Perm_q_CT2 = qvals_CT2;

% Count discoveries
R_CT2 = sum(qvals_CT2(validIdx_perm_CT2) < 0.03);

% Calculate FCR-adjusted alpha
if R_CT2 == 0
    alpha_adj_CT2 = 0.05 / m_CT2;
    fprintf('CT2 FCR: No discoveries. Using Bonferroni.\n');
else
    alpha_adj_CT2 = 0.05 * (R_CT2 / m_CT2);
    fprintf('CT2 FCR: %d discoveries. Adjusted alpha = %.5f.\n', R_CT2, alpha_adj_CT2);
end

% Adjust bootstrap BCa CIs using FCR
Boot_CIlo_bca_fcr_CT2 = nan(nAreas, 1);
Boot_CIhi_bca_fcr_CT2 = nan(nAreas, 1);
Boot_CIlo_bca_bonf_CT2 = nan(nAreas, 1);
Boot_CIhi_bca_bonf_CT2 = nan(nAreas, 1);

for i = 1:nAreas
    if ~isempty(Boot_distributions_CT2{i}) && ~isnan(Boot_z0_CT2(i))
        z0 = Boot_z0_CT2(i);
        a = Boot_a_CT2(i);
        
        % FCR adjustment
        z_alpha_lower_fcr = norminv(alpha_adj_CT2 / 2);
        z_alpha_upper_fcr = norminv(1 - alpha_adj_CT2 / 2);
        
        alpha_lower_adj_fcr = normcdf(z0 + (z0 + z_alpha_lower_fcr) / (1 - a * (z0 + z_alpha_lower_fcr)));
        alpha_upper_adj_fcr = normcdf(z0 + (z0 + z_alpha_upper_fcr) / (1 - a * (z0 + z_alpha_upper_fcr)));
        
        alpha_lower_adj_fcr = max(0.001, min(0.999, alpha_lower_adj_fcr));
        alpha_upper_adj_fcr = max(0.001, min(0.999, alpha_upper_adj_fcr));
        
        ci_fcr = prctile(Boot_distributions_CT2{i}, [alpha_lower_adj_fcr * 100, alpha_upper_adj_fcr * 100]);
        Boot_CIlo_bca_fcr_CT2(i) = ci_fcr(1);
        Boot_CIhi_bca_fcr_CT2(i) = ci_fcr(2);
        
        % Bonferroni adjustment
        z_alpha_lower_bonf = norminv(alpha_bonf / 2);
        z_alpha_upper_bonf = norminv(1 - alpha_bonf / 2);
        
        alpha_lower_adj_bonf = normcdf(z0 + (z0 + z_alpha_lower_bonf) / (1 - a * (z0 + z_alpha_lower_bonf)));
        alpha_upper_adj_bonf = normcdf(z0 + (z0 + z_alpha_upper_bonf) / (1 - a * (z0 + z_alpha_upper_bonf)));
        
        alpha_lower_adj_bonf = max(0.001, min(0.999, alpha_lower_adj_bonf));
        alpha_upper_adj_bonf = max(0.001, min(0.999, alpha_upper_adj_bonf));
        
        ci_bonf = prctile(Boot_distributions_CT2{i}, [alpha_lower_adj_bonf * 100, alpha_upper_adj_bonf * 100]);
        Boot_CIlo_bca_bonf_CT2(i) = ci_bonf(1);
        Boot_CIhi_bca_bonf_CT2(i) = ci_bonf(2);
    end
end

% --- Interaction FCR ---
validIdx_perm_int = ~isnan(Perm_int_p);
m_int = sum(validIdx_perm_int);

% Calculate q-values based on permutation p-values
qvals_int = nan(nAreas, 1);
if m_int > 0
    qvals_int(validIdx_perm_int) = mafdr(Perm_int_p(validIdx_perm_int), 'BHFDR', true);
end
Perm_int_q = qvals_int;

% Count discoveries
R_int = sum(qvals_int(validIdx_perm_int) < 0.045);

% Calculate FCR-adjusted alpha
if R_int == 0
    alpha_adj_int = 0.05 / m_int;
    fprintf('Interaction FCR: No discoveries. Using Bonferroni.\n');
else
    alpha_adj_int = 0.05 * (R_int / m_int);
    fprintf('Interaction FCR: %d discoveries. Adjusted alpha = %.5f.\n', R_int, alpha_adj_int);
end

% Adjust bootstrap BCa CIs using FCR
Boot_int_CIlo_bca_fcr = nan(nAreas, 1);
Boot_int_CIhi_bca_fcr = nan(nAreas, 1);
Boot_int_CIlo_bca_bonf = nan(nAreas, 1);
Boot_int_CIhi_bca_bonf = nan(nAreas, 1);

for i = 1:nAreas
    if ~isempty(Boot_int_distributions{i}) && ~isnan(Boot_int_z0(i))
        z0 = Boot_int_z0(i);
        a = Boot_int_a(i);
        
        % FCR adjustment
        z_alpha_lower_fcr = norminv(alpha_adj_int / 2);
        z_alpha_upper_fcr = norminv(1 - alpha_adj_int / 2);
        
        alpha_lower_adj_fcr = normcdf(z0 + (z0 + z_alpha_lower_fcr) / (1 - a * (z0 + z_alpha_lower_fcr)));
        alpha_upper_adj_fcr = normcdf(z0 + (z0 + z_alpha_upper_fcr) / (1 - a * (z0 + z_alpha_upper_fcr)));
        
        alpha_lower_adj_fcr = max(0.001, min(0.999, alpha_lower_adj_fcr));
        alpha_upper_adj_fcr = max(0.001, min(0.999, alpha_upper_adj_fcr));
        
        ci_fcr = prctile(Boot_int_distributions{i}, [alpha_lower_adj_fcr * 100, alpha_upper_adj_fcr * 100]);
        Boot_int_CIlo_bca_fcr(i) = ci_fcr(1);
        Boot_int_CIhi_bca_fcr(i) = ci_fcr(2);
        
        % Bonferroni adjustment
        z_alpha_lower_bonf = norminv(alpha_bonf / 2);
        z_alpha_upper_bonf = norminv(1 - alpha_bonf / 2);
        
        alpha_lower_adj_bonf = normcdf(z0 + (z0 + z_alpha_lower_bonf) / (1 - a * (z0 + z_alpha_lower_bonf)));
        alpha_upper_adj_bonf = normcdf(z0 + (z0 + z_alpha_upper_bonf) / (1 - a * (z0 + z_alpha_upper_bonf)));
        
        alpha_lower_adj_bonf = max(0.001, min(0.999, alpha_lower_adj_bonf));
        alpha_upper_adj_bonf = max(0.001, min(0.999, alpha_upper_adj_bonf));
        
        ci_bonf = prctile(Boot_int_distributions{i}, [alpha_lower_adj_bonf * 100, alpha_upper_adj_bonf * 100]);
        Boot_int_CIlo_bca_bonf(i) = ci_bonf(1);
        Boot_int_CIhi_bca_bonf(i) = ci_bonf(2);
    end
end

% ==================== Summary ====================
fprintf('\n=======================================================\n');
fprintf('BOOTSTRAP + PERMUTATION SUMMARY\n');
fprintf('=======================================================\n');
fprintf('CT1:\n');
fprintf('  Permutation FDR discoveries: %d\n', R_CT1);
fprintf('  Bootstrap BCa (unadj) CI excl. 0: %d\n', sum((Boot_CIlo_bca_CT1 > 0) | (Boot_CIhi_bca_CT1 < 0)));
fprintf('  Bootstrap BCa (FCR) CI excl. 0: %d\n', sum((Boot_CIlo_bca_fcr_CT1 > 0) | (Boot_CIhi_bca_fcr_CT1 < 0)));
fprintf('  Bootstrap BCa (Bonf) CI excl. 0: %d\n', sum((Boot_CIlo_bca_bonf_CT1 > 0) | (Boot_CIhi_bca_bonf_CT1 < 0)));
fprintf('CT2:\n');
fprintf('  Permutation FDR discoveries: %d\n', R_CT2);
fprintf('  Bootstrap BCa (unadj) CI excl. 0: %d\n', sum((Boot_CIlo_bca_CT2 > 0) | (Boot_CIhi_bca_CT2 < 0)));
fprintf('  Bootstrap BCa (FCR) CI excl. 0: %d\n', sum((Boot_CIlo_bca_fcr_CT2 > 0) | (Boot_CIhi_bca_fcr_CT2 < 0)));
fprintf('  Bootstrap BCa (Bonf) CI excl. 0: %d\n', sum((Boot_CIlo_bca_bonf_CT2 > 0) | (Boot_CIhi_bca_bonf_CT2 < 0)));
fprintf('Interaction:\n');
fprintf('  Permutation FDR discoveries: %d\n', R_int);
fprintf('  Bootstrap BCa (unadj) CI excl. 0: %d\n', sum((Boot_int_CIlo_bca > 0) | (Boot_int_CIhi_bca < 0)));
fprintf('  Bootstrap BCa (FCR) CI excl. 0: %d\n', sum((Boot_int_CIlo_bca_fcr > 0) | (Boot_int_CIhi_bca_fcr < 0)));
fprintf('  Bootstrap BCa (Bonf) CI excl. 0: %d\n', sum((Boot_int_CIlo_bca_bonf > 0) | (Boot_int_CIhi_bca_bonf < 0)));
fprintf('=======================================================\n');
% % -------------------- Build Results table --------------------
% Results = table( ...
%     areaNames(:), ...
%     meanFracU_CT1, meanFracL_CT1, log2FC_CT1, ...
%     GLMM_eff_CT1, GLMM_se_CT1, GLMM_CIlo_CT1, GLMM_CIhi_CT1, ...
%     GLMM_CIlo_CT1_bonf, GLMM_CIhi_CT1_bonf, GLMM_CIlo_CT1_fcr, GLMM_CIhi_CT1_fcr, ...
%     GLMM_p_CT1, q_GLMM_CT1, ranksum_p_CT1, ...
%     meanFracU_CT2, meanFracL_CT2, log2FC_CT2, ...
%     GLMM_eff_CT2, GLMM_se_CT2, GLMM_CIlo_CT2, GLMM_CIhi_CT2, ...
%     GLMM_CIlo_CT2_bonf, GLMM_CIhi_CT2_bonf, GLMM_CIlo_CT2_fcr, GLMM_CIhi_CT2_fcr, ...
%     GLMM_p_CT2, q_GLMM_CT2, ranksum_p_CT2, ...
%     GLMM_int_eff, GLMM_int_se, GLMM_int_CIlo, GLMM_int_CIhi, ...
%     GLMM_int_CIlo_bonf, GLMM_int_CIhi_bonf, GLMM_int_CIlo_fcr, GLMM_int_CIhi_fcr, ...
%     GLMM_int_p, GLMM_int_q, ...
%     Boot_int_eff, Boot_int_CIlo_bca, Boot_int_CIhi_bca, ...
%     Boot_int_CIlo_bca_fcr, Boot_int_CIhi_bca_fcr, ...
%     Boot_int_CIlo_bca_bonf, Boot_int_CIhi_bca_bonf, ...
%     Boot_int_p, Boot_int_q, ...
%     'VariableNames', {'Area', ...
%     'meanFracU_CT1','meanFracL_CT1','log2FC_CT1', ...
%     'GLMM_logOdds_U_vs_L_CT1','GLMM_SE_CT1','GLMM_CIlo_CT1','GLMM_CIhi_CT1', ...
%     'GLMM_CIlo_CT1_Bonf','GLMM_CIhi_CT1_Bonf','GLMM_CIlo_CT1_FCR','GLMM_CIhi_CT1_FCR', ...
%     'GLMM_p_CT1','GLMM_q_CT1','ranksum_p_CT1', ...
%     'meanFracU_CT2','meanFracL_CT2','log2FC_CT2', ...
%     'GLMM_logOdds_U_vs_L_CT2','GLMM_SE_CT2','GLMM_CIlo_CT2','GLMM_CIhi_CT2', ...
%     'GLMM_CIlo_CT2_Bonf','GLMM_CIhi_CT2_Bonf','GLMM_CIlo_CT2_FCR','GLMM_CIhi_CT2_FCR', ...
%     'GLMM_p_CT2','GLMM_q_CT2','ranksum_p_CT2', ...
%     'Interaction_logOddsDiff','Interaction_SE','Interaction_CIlo','Interaction_CIhi', ...
%     'Interaction_CIlo_Bonf','Interaction_CIhi_Bonf','Interaction_CIlo_FCR','Interaction_CIhi_FCR', ...
%     'Interaction_p','Interaction_q',...
%     'Boot_Interaction_logOddsDiff', 'Boot_Interaction_CIlo_BCa', 'Boot_Interaction_CIhi_BCa', ...
%     'Boot_Interaction_CIlo_BCa_FCR', 'Boot_Interaction_CIhi_BCa_FCR', ...
%     'Boot_Interaction_CIlo_BCa_Bonf', 'Boot_Interaction_CIhi_BCa_Bonf', ...
%     'Boot_Interaction_p', 'Boot_Interaction_q'});


% -------------------- Build Results table --------------------
Results = table(areaNames(:), ...
    meanFracU_CT1, meanFracL_CT1, log2FC_CT1, ...
    Boot_eff_CT1, Boot_CIlo_bca_CT1, Boot_CIhi_bca_CT1, ...
    Boot_CIlo_bca_fcr_CT1, Boot_CIhi_bca_fcr_CT1, ...
    Boot_CIlo_bca_bonf_CT1, Boot_CIhi_bca_bonf_CT1, ...
    Perm_p_CT1, Perm_q_CT1, ...
    meanFracU_CT2, meanFracL_CT2, log2FC_CT2, ...
    Boot_eff_CT2, Boot_CIlo_bca_CT2, Boot_CIhi_bca_CT2, ...
    Boot_CIlo_bca_fcr_CT2, Boot_CIhi_bca_fcr_CT2, ...
    Boot_CIlo_bca_bonf_CT2, Boot_CIhi_bca_bonf_CT2, ...
    Perm_p_CT2, Perm_q_CT2, ...
    Boot_int_eff, Boot_int_CIlo_bca, Boot_int_CIhi_bca, ...
    Boot_int_CIlo_bca_fcr, Boot_int_CIhi_bca_fcr, ...
    Boot_int_CIlo_bca_bonf, Boot_int_CIhi_bca_bonf, ...
    Perm_int_p, Perm_int_q, ...
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
    'Perm_p_CT2', 'Perm_q_CT2',...
    'Boot_Interaction', 'Boot_Int_CIlo_BCa', 'Boot_Int_CIhi_BCa', ...
    'Boot_Int_CIlo_BCa_FCR', 'Boot_Int_CIhi_BCa_FCR', ...
    'Boot_Int_CIlo_BCa_Bonf', 'Boot_Int_CIhi_BCa_Bonf', ...
    'Perm_Int_p', 'Perm_Int_q'});


% % Overwrite unadjusted CIs with FCR (or comment out to keep unadjusted)
% Results.GLMM_CIlo_CT1 = GLMM_CIlo_CT1_fcr;
% Results.GLMM_CIhi_CT1 = GLMM_CIhi_CT1_fcr;
% Results.GLMM_CIlo_CT2 = GLMM_CIlo_CT2_fcr;
% Results.GLMM_CIhi_CT2 = GLMM_CIhi_CT2_fcr;
% Results.Interaction_CIlo = GLMM_int_CIlo_fcr;
% Results.Interaction_CIhi = GLMM_int_CIhi_fcr;

% % Overwrite unadjusted CIs with Bonferroni (or comment out to keep unadjusted)
% Results.GLMM_CIlo_CT1 = GLMM_CIlo_CT1_bonf;
% Results.GLMM_CIhi_CT1 = GLMM_CIhi_CT1_bonf;
% Results.GLMM_CIlo_CT2 = GLMM_CIlo_CT2_bonf;
% Results.GLMM_CIhi_CT2 = GLMM_CIhi_CT2_bonf;
% Results.Interaction_CIlo = GLMM_int_CIlo_bonf;
% Results.Interaction_CIhi = GLMM_int_CIhi_bonf;

% -------------------- PLOTTING (3 panels) --------------------
if doPlot
    topN = 6;
    nAreasTotal = numel(Results.Area);

    % Abundance fractions
    abund_frac_ct1 = 0.5*(Results.meanFracU_CT1 + Results.meanFracL_CT1);
    abund_frac_ct2 = 0.5*(Results.meanFracU_CT2 + Results.meanFracL_CT2);
    abund_frac_tot  = 0.25*((Results.meanFracU_CT1 + Results.meanFracL_CT1) + (Results.meanFracU_CT2 + Results.meanFracL_CT2));

    v1 = 100 * abund_frac_ct1; v1(isnan(v1)) = 0;
    v2 = 100 * abund_frac_ct2; v2(isnan(v2)) = 0;
    vtot = 100 * abund_frac_tot; vtot(isnan(vtot)) = 0;

    % Build INTERACTION ordering (CHANGED FROM CT1)
    effInt_all = Results.Boot_Interaction;
    effTmp = effInt_all;
    effTmp(isnan(effTmp)) = -Inf;
    [~, orderInt] = sort(effTmp, 'descend');

    % Create figure
    figA = figure('Units','centimeters','Position',[2 2 34 18],'Color','w');
    tA = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    % Size scaling
    sMin = 8; sMax = 260;
    if max(v1) == min(v1)
        sVals1 = repmat(mean([sMin sMax]), nAreasTotal, 1);
    else
        sVals1 = sMin + (v1 - min(v1))./(max(v1)-min(v1))*(sMax-sMin);
    end

    % CT1 panel
    ax1 = nexttile(tA,1);
    draw_forest_aligned(ax1, Results.Boot_logOdds_CT1, ...
        Results.Boot_CIlo_BCa_FCR_CT1, Results.Boot_CIhi_BCa_FCR_CT1, ...
        Results.Area, Results.Perm_q_CT1, ...
        ColorIntPos, ColorIntPos, ColorIntPos, sVals1, 'CT1: log-odds (U vs L)', orderInt, topN); hold on;

    % draw_forest_aligned(ax1, Results.Boot_logOdds_CT1, ...
    %     Results.Boot_CIlo_BCa_FCR_CT1, Results.Boot_CIhi_BCa_FCR_CT1, ...
    %     Results.Area, Results.Perm_q_CT1, ...
    %     ColorU1, ColorL1, ColorNonSig, sVals1, 'CT1: log-odds (U vs L)', orderCT1, topN); hold on;

    if isempty(v1) || all(isnan(v1))
        vals = [10, 5, 1];
    else
        posOnly = v1(v1 > 0);
        if isempty(posOnly)
            vals = [max(v1), mean(v1), min(v1)];
        else
            vals = [max(v1), mean(posOnly), min(posOnly)];
        end
    end

    normL = (vals - min(v1)) / (max(v1)-min(v1));
    sThree = sMin + normL * (sMax - sMin);

    yspot = [27, 28, 29];
    xlims=xlim;
    for ii = 1:3
        scatter(ax1, xlims(2)/2, yspot(ii), sThree(ii), 'o', 'MarkerFaceColor', ColorNonSig, 'MarkerEdgeColor', ColorNonSig);
        text(ax1, xlims(2)/2+0.5, yspot(ii), sprintf('%.1f%%', vals(ii)), 'FontSize', 11, 'VerticalAlignment', 'middle');
    end
    hold off;

    %size scaling

    sMin = 8; sMax = 260;
    if max(v2) == min(v2)
        sVals2 = repmat(mean([sMin sMax]), nAreasTotal, 1);
    else
        sVals2 = sMin + (v2 - min(v2))./(max(v2)-min(v2))*(sMax-sMin);
    end

    % CT2 panel
    ax2 = nexttile(tA,2);
    draw_forest_aligned(ax2, Results.Boot_logOdds_CT2, ...
        Results.Boot_CIlo_BCa_FCR_CT2, Results.Boot_CIhi_BCa_FCR_CT2, ...
        Results.Area, Results.Perm_q_CT2, ...
        ColorIntNeg, ColorIntNeg, ColorIntNeg, sVals2, 'CT2: log-odds (U vs L)', orderInt, topN); hold on;

    % draw_forest_aligned(ax2, Results.Boot_logOdds_CT2, ...
    %     Results.Boot_CIlo_BCa_FCR_CT2, Results.Boot_CIhi_BCa_FCR_CT2, ...
    %     Results.Area, Results.Perm_q_CT2, ...
    %     ColorU2, ColorL2, ColorNonSig, sVals2, 'CT2: log-odds (U vs L)', orderCT1, topN); hold on;



    if isempty(v2) || all(isnan(v2))
        vals = [10, 5, 1];
    else
        posOnly = v2(v2 > 0);
        if isempty(posOnly)
            vals = [max(v2), mean(v2), min(v2)];
        else
            vals = [max(v2), mean(posOnly), min(posOnly)];
        end
    end

    normL = (vals - min(v2)) / (max(v2)-min(v2));
    sThree = sMin + normL * (sMax - sMin);

    yspot = [27, 28, 29];
    xlims=xlim;
    for ii = 1:3
        scatter(ax2, xlims(2)/2, yspot(ii), sThree(ii), 'o', 'MarkerFaceColor', ColorNonSig, 'MarkerEdgeColor', ColorNonSig);
        text(ax2, xlims(2)/2+0.7, yspot(ii), sprintf('%.1f%%', vals(ii)), 'FontSize', 11, 'VerticalAlignment', 'middle');
    end
    hold off;

    %size scaling
    sMin = 8; sMax = 260;
    if max(vtot) == min(vtot)
        sValstot = repmat(mean([sMin sMax]), nAreasTotal, 1);
    else
        sValstot = sMin + (vtot - min(vtot))./(max(vtot)-min(vtot))*(sMax-sMin);
    end

    % Interaction panel (Bootstrap BCa with FCR based on permutation)
    ax3 = nexttile(tA,3);
    draw_forest_aligned(ax3, Results.Boot_Interaction, ...
        Results.Boot_Int_CIlo_BCa_FCR, Results.Boot_Int_CIhi_BCa_FCR, ...
        Results.Area, Results.Perm_Int_q, ...
        [0.4 0.4 0.4], [0.4 0.4 0.4], [0.4 0.4 0.4], sValstot, 'CT1-CT2 interaction', orderInt, topN); hold on;

    % draw_forest_aligned(ax3, Results.Boot_Interaction, ...
    %     Results.Boot_Int_CIlo_BCa_FCR, Results.Boot_Int_CIhi_BCa_FCR, ...
    %     Results.Area, Results.Perm_Int_q, ...
    %     ColorIntPos, ColorIntNeg, ColorNonSig, sValstot, 'CT1-CT2 interaction', orderCT1, topN); hold on;

    if isempty(vtot) || all(isnan(vtot))
        vals = [10, 5, 1];
    else
        posOnly = vtot(vtot > 0);
        if isempty(posOnly)
            vals = [max(vtot), mean(vtot), min(vtot)];
        else
            vals = [max(vtot), mean(posOnly), min(posOnly)];
        end
    end

    normL = (vals - min(vtot)) / (max(vtot)-min(vtot));
    sThree = sMin + normL * (sMax - sMin);

    yspot = [27, 28, 29];
    xlims=xlim;
    for ii = 1:3
        scatter(ax3, xlims(2)/2, yspot(ii), sThree(ii), 'o', 'MarkerFaceColor', ColorNonSig, 'MarkerEdgeColor', ColorNonSig);
        text(ax3, xlims(2)/2+0.3, yspot(ii), sprintf('%.1f%%', vals(ii)), 'FontSize', 11, 'VerticalAlignment', 'middle');
    end

    hold off;

    % % Volcano Plot for Interaction (Effect Size vs -log10(p-value))
    % figure('Units','centimeters','Position',[1 15 20.4047 8],'Color','w');
    % 
    % valid = ~isnan(Perm_int_p) & ~isnan(Boot_int_eff);
    % 
    % if any(valid)
    %     effs = Boot_int_eff(valid);
    %     pvals = Perm_int_p(valid);
    %     qvals = Perm_int_q(valid);
    %     names = areaNames(valid);
    % 
    %     % Calculate -log10(p)
    %     neg_log_p = -log10(pvals);
    % 
    %     % Thresholds
    %     fdr_threshold = 0.056;
    % 
    %     hold on;
    % 
    %     % Plot each point
    %     for k = 1:numel(effs)
    %         e = effs(k);
    %         nlp = neg_log_p(k);
    %         q = qvals(k);
    % 
    % 
    %         % Determine color and size
    %         if q < fdr_threshold
    %             if e > 0
    %                 col = ColorIntPos;  % CT1 > CT2
    %             else
    %                 col = ColorIntNeg;  % CT2 > CT1
    %             end
    %             sz = 80;
    %         else
    %             col = ColorNonSig;
    %             sz = 40;
    %         end
    % 
    %         scatter(e, nlp, sz, 'o', 'MarkerFaceColor', col, ...
    %             'MarkerEdgeColor', col*0.6, 'LineWidth', 0.5);
    % 
    %         % Only label significant points or top effects
    %         if q < fdr_threshold
    %             text(e+0.05, nlp, names{k}, ...
    %                 'FontSize', 9, 'HorizontalAlignment', 'left', ...
    %                 'VerticalAlignment', 'middle');
    %         end
    %     end
    % 
    %     % Add threshold lines
    %     yline(-log10(0.05), '--', 'p = 0.05', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    %     xline(0, '-k', 'LineWidth', 1);
    % 
    %     % Labels
    %     xlabel('Interaction: log-odds difference (CT1 - CT2)', 'FontSize', 11);
    %     ylabel('-log_{10}(p-value)', 'FontSize', 11);
    %     title('Interaction Permutation Test: Volcano Plot', 'FontSize', 12);
    % 
    %     grid on;
    %     ax = gca;
    %     ax.GridColor = [0.2 0.2 0.2];
    %     ax.GridAlpha = 0.1;
    %     box off;
    % 
    %     hold off;
    % else
    %     text(0.5, 0.5, 'No valid interaction results to plot', ...
    %         'FontSize', 11, 'HorizontalAlignment', 'center');
    % end
    % 

    %--------------3 volcano plot panels ----------------

    figure('Units','centimeters','Position',[1 5 45 10],'Color','w');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    fdr_threshold = 0.05;

    % ==================== CT1 Volcano ====================
    ax1 = nexttile(1);
    valid_CT1 = ~isnan(Perm_p_CT1) & ~isnan(Boot_eff_CT1);

    if any(valid_CT1)
        effs = Boot_eff_CT1(valid_CT1);
        pvals = Perm_p_CT1(valid_CT1);
        qvals = Perm_q_CT1(valid_CT1);
        names = areaNames(valid_CT1);
        neg_log_p = -log10(pvals);

        hold on;
        for k = 1:numel(effs)
            e = effs(k);
            nlp = neg_log_p(k);
            q = qvals(k);
            p=pvals(k);

            if q < fdr_threshold
                if e > 0
                    col = ColorU1;
                else
                    col = ColorL1;
                end
                sz = 20;
            else
                col = ColorNonSig;
                sz = 20;
            end

            if p<0.052
                text(e+0.05, nlp, names{k}, 'FontSize', 8, 'HorizontalAlignment', 'left');
            end

            scatter(e, nlp, sz, 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', col*0.6, 'LineWidth', 0.5);
        end

        yline(-log10(0.05), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        xline(0, '-k', 'LineWidth', 1);
        xlabel('Log-odds ratio (U vs L)', 'FontSize', 10);
        ylabel('-log_{10}(p-value)', 'FontSize', 10);
        title('CT1: Volcano Plot', 'FontSize', 11);
        grid on; ax1.GridAlpha = 0.1; box off;
        hold off;
    end

    % ==================== CT2 Volcano ====================
    ax2 = nexttile(2);
    valid_CT2 = ~isnan(Perm_p_CT2) & ~isnan(Boot_eff_CT2);

    if any(valid_CT2)
        effs = Boot_eff_CT2(valid_CT2);
        pvals = Perm_p_CT2(valid_CT2);
        qvals = Perm_q_CT2(valid_CT2);
        names = areaNames(valid_CT2);
        neg_log_p = -log10(pvals);

        hold on;
        for k = 1:numel(effs)
            e = effs(k);
            nlp = neg_log_p(k);
            q = qvals(k);
            p = pvals(k);

            if q < 0.04 %fdr_threshold
                if e > 0
                    col = ColorU2;
                else
                    col = ColorL2;
                end
                sz = 20;
                
            else
                col = ColorNonSig;
                sz = 20;
            end
            
            if p<0.05
                text(e+0.05, nlp, names{k}, 'FontSize', 8, 'HorizontalAlignment', 'left');
            end

            scatter(e, nlp, sz, 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', col*0.6, 'LineWidth', 0.5);
        end

        yline(-log10(0.05), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        xline(0, '-k', 'LineWidth', 1);
        xlabel('Log-odds ratio (U vs L)', 'FontSize', 10);
        ylabel('-log_{10}(p-value)', 'FontSize', 10);
        title('CT2: Volcano Plot', 'FontSize', 11);
        grid on; ax2.GridAlpha = 0.1; box off;
        hold off;
    end

    % ==================== Interaction Volcano ====================
    ax3 = nexttile(3);
    valid_int = ~isnan(Perm_int_p) & ~isnan(Boot_int_eff);

    if any(valid_int)
        effs = Boot_int_eff(valid_int);
        pvals = Perm_int_p(valid_int);
        qvals = Perm_int_q(valid_int);
        names = areaNames(valid_int);
        neg_log_p = -log10(pvals);

        hold on;
        for k = 1:numel(effs)
            e = effs(k);
            nlp = neg_log_p(k);
            q = qvals(k);
            p = pvals(k);

            if p < 0.05 %fdr_threshold
                if e > 0 && q<fdr_threshold
                    col = ColorIntPos;
                elseif e<=0 && q<fdr_threshold
                    col = ColorIntNeg;
                end
                sz = 20;
                text(e+0.05, nlp, names{k}, 'FontSize', 8, 'HorizontalAlignment', 'left');
            else
                col = ColorNonSig;
                sz = 20;
            end

            scatter(e, nlp, sz, 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', col*0.6, 'LineWidth', 0.5);
        end

        yline(-log10(0.05), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        xline(0, '-k', 'LineWidth', 1);
        xlabel('Interaction (CT1 - CT2)', 'FontSize', 10);
        ylabel('-log_{10}(p-value)', 'FontSize', 10);
        title('Interaction: Volcano Plot', 'FontSize', 11);
        grid on; ax3.GridAlpha = 0.1; box off;
        hold off;
    end
    % topN = 6;
    % nAreasTotal = numel(Results.Area);
    %
    % % Abundance fractions
    % abund_frac_ct1 = 0.5*(Results.meanFracU_CT1 + Results.meanFracL_CT1);
    % abund_frac_ct2 = 0.5*(Results.meanFracU_CT2 + Results.meanFracL_CT2);
    % abund_frac_tot  = 0.25*((Results.meanFracU_CT1 + Results.meanFracL_CT1) + (Results.meanFracU_CT2 + Results.meanFracL_CT2));
    % 
    % v1 = 100 * abund_frac_ct1; v1(isnan(v1)) = 0;
    % v2 = 100 * abund_frac_ct2; v2(isnan(v2)) = 0;
    % vtot = 100 * abund_frac_tot; vtot(isnan(vtot)) = 0;
    % 
    % % Build CT1 ordering
    % effCT1_all = Results.GLMM_logOdds_U_vs_L_CT1;
    % effTmp = effCT1_all;
    % effTmp(isnan(effTmp)) = -Inf;
    % [~, orderCT1] = sort(effTmp, 'descend');
    % 
    % % Create figure
    % figA = figure('Units','centimeters','Position',[2 2 34 18],'Color','w');
    % tA = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    % 
    % sMin = 8; sMax = 260;
    % if max(v1) == min(v1)
    %     sVals1 = repmat(mean([sMin sMax]), nAreasTotal, 1);
    % else
    %     sVals1 = sMin + (v1 - min(v1))./(max(v1)-min(v1))*(sMax-sMin);
    % end
    % 
    % % CT1 panel
    % ax1 = nexttile(tA,1);
    % draw_forest_aligned(ax1, Results.GLMM_logOdds_U_vs_L_CT1, Results.GLMM_CIlo_CT1, Results.GLMM_CIhi_CT1, Results.Area, ...
    %     Results.GLMM_q_CT1, ColorU1, ColorL1, ColorNonSig, sVals1, 'CT1: log-odds (U vs L)', orderCT1, topN); hold on;
    % 
    % if isempty(v1) || all(isnan(v1))
    %     vals = [10, 5, 1];
    % else
    %     posOnly = v1(v1 > 0);
    %     if isempty(posOnly)
    %         vals = [max(v1), mean(v1), min(v1)];
    %     else
    %         vals = [max(v1), mean(posOnly), min(posOnly)];
    %     end
    % end
    % 
    % normL = (vals - min(v1)) / (max(v1)-min(v1));
    % sThree = sMin + normL * (sMax - sMin);
    % 
    % yspot = [27, 28, 29];
    % xlims=xlim;
    % for ii = 1:3
    %     scatter(ax1, xlims(2)/2, yspot(ii), sThree(ii), 'o', 'MarkerFaceColor', ColorNonSig, 'MarkerEdgeColor', ColorNonSig);
    %     text(ax1, xlims(2)/2+0.5, yspot(ii), sprintf('%.1f%%', vals(ii)), 'FontSize', 11, 'VerticalAlignment', 'middle');
    % end
    % hold off;
    % 
    % sMin = 8; sMax = 260;
    % if max(v2) == min(v2)
    %     sVals2 = repmat(mean([sMin sMax]), nAreasTotal, 1);
    % else
    %     sVals2 = sMin + (v2 - min(v2))./(max(v2)-min(v2))*(sMax-sMin);
    % end
    % 
    % % CT2 panel
    % ax2 = nexttile(tA,2);
    % draw_forest_aligned(ax2, Results.GLMM_logOdds_U_vs_L_CT2, Results.GLMM_CIlo_CT2, Results.GLMM_CIhi_CT2, Results.Area, ...
    %     Results.GLMM_q_CT2, ColorU2, ColorL2, ColorNonSig, sVals2, 'CT2: log-odds (U vs L)', orderCT1, topN); hold on;
    % 
    % if isempty(v2) || all(isnan(v2))
    %     vals = [10, 5, 1];
    % else
    %     posOnly = v2(v2 > 0);
    %     if isempty(posOnly)
    %         vals = [max(v2), mean(v2), min(v2)];
    %     else
    %         vals = [max(v2), mean(posOnly), min(posOnly)];
    %     end
    % end
    % 
    % normL = (vals - min(v2)) / (max(v2)-min(v2));
    % sThree = sMin + normL * (sMax - sMin);
    % 
    % yspot = [27, 28, 29];
    % xlims=xlim;
    % for ii = 1:3
    %     scatter(ax2, xlims(2)/2, yspot(ii), sThree(ii), 'o', 'MarkerFaceColor', ColorNonSig, 'MarkerEdgeColor', ColorNonSig);
    %     text(ax2, xlims(2)/2+0.7, yspot(ii), sprintf('%.1f%%', vals(ii)), 'FontSize', 11, 'VerticalAlignment', 'middle');
    % end
    % hold off;
    % 
    % sMin = 8; sMax = 260;
    % if max(vtot) == min(vtot)
    %     sValstot = repmat(mean([sMin sMax]), nAreasTotal, 1);
    % else
    %     sValstot = sMin + (vtot - min(vtot))./(max(vtot)-min(vtot))*(sMax-sMin);
    % end
    % 
    % % Interaction panel
    % % ax3 = nexttile(tA,3);
    % % draw_forest_aligned(ax3, Results.Interaction_logOddsDiff, Results.Interaction_CIlo, Results.Interaction_CIhi, Results.Area, ...
    % %     Results.Interaction_q, ColorIntPos, ColorIntNeg, ColorNonSig, sValstot, 'CT1-CT2 interaction: log-odds difference', orderCT1, topN); hold on;
    % 
    % ax3 = nexttile(tA,3);
    % draw_forest_aligned(ax3, Results.Boot_Interaction_logOddsDiff, ...
    %     Results.Boot_Interaction_CIlo_BCa_FCR, Results.Boot_Interaction_CIhi_BCa_FCR, ...
    %     Results.Area, Results.Boot_Interaction_q, ...
    %     ColorIntPos, ColorIntNeg, ColorNonSig, sValstot, ...
    % 'CT1-CT2 interaction: log-odds difference (Bootstrap BCa)', orderCT1, topN);
    % hold on;
    % 
    % if isempty(vtot) || all(isnan(vtot))
    %     vals = [10, 5, 1];
    % else
    %     posOnly = vtot(vtot > 0);
    %     if isempty(posOnly)
    %         vals = [max(vtot), mean(vtot), min(vtot)];
    %     else
    %         vals = [max(vtot), mean(posOnly), min(posOnly)];
    %     end
    % end
    % 
    % normL = (vals - min(vtot)) / (max(vtot)-min(vtot));
    % sThree = sMin + normL * (sMax - sMin);
    % 
    % yspot = [27, 28, 29];
    % xlims=xlim;
    % for ii = 1:3
    %     scatter(ax3, xlims(2)/2, yspot(ii), sThree(ii), 'o', 'MarkerFaceColor', ColorNonSig, 'MarkerEdgeColor', ColorNonSig);
    %     text(ax3, xlims(2)/2+0.3, yspot(ii), sprintf('%.1f%%', vals(ii)), 'FontSize', 11, 'VerticalAlignment', 'middle');
    % end
    % 
    % hold off;

    % % ---------------- FIG B: CT1 vs CT2 log-odds correlation with BCa CIs and names----------------
    % 
    % x_est = Boot_eff_CT1;
    % y_est = Boot_eff_CT2;
    % 
    % % Use FCR-corrected BCa CIs from the main analysis as error bars
    % x_err_lo = Boot_eff_CT1 - Boot_CIlo_bca_fcr_CT1;
    % x_err_hi = Boot_CIhi_bca_fcr_CT1 - Boot_eff_CT1;
    % y_err_lo = Boot_eff_CT2 - Boot_CIlo_bca_fcr_CT2;
    % y_err_hi = Boot_CIhi_bca_fcr_CT2 - Boot_eff_CT2;
    % 
    % validPair = isfinite(x_est) & isfinite(y_est);
    % 
    % if sum(validPair) > 2
    %     xv = x_est(validPair);
    %     yv = y_est(validPair);
    %     xel = x_err_lo(validPair);
    %     xeh = x_err_hi(validPair);
    %     yel = y_err_lo(validPair);
    %     yeh = y_err_hi(validPair);
    %     areasPair = areaNames(validPair);
    % 
    %     [r_obs, p_obs] = corr(xv, yv, 'type', 'Pearson');
    % 
    %     nboot_corr = 10000;
    %     boot_r_corr = zeros(nboot_corr, 1);
    % 
    %     x_fit_corr = linspace(min(xv), max(xv), 100);
    %     y_boot_corr = zeros(nboot_corr, length(x_fit_corr));
    % 
    %     rng(42);
    %     for b = 1:nboot_corr
    %         idx_U1 = randi(nU1, nU1, 1);
    %         idx_L1 = randi(nL1, nL1, 1);
    %         idx_U2 = randi(nU2, nU2, 1);
    %         idx_L2 = randi(nL2, nL2, 1);
    % 
    %         boot_logOR_CT1 = nan(nAreas, 1);
    %         boot_logOR_CT2 = nan(nAreas, 1);
    % 
    %         for i = 1:nAreas
    %             succ_U1_b = nu1(idx_U1, i);
    %             succ_L1_b = nl1(idx_L1, i);
    %             trials_U1_b = N_U1(idx_U1);
    %             trials_L1_b = N_L1(idx_L1);
    % 
    %             if ~any(isnan(succ_U1_b)) && ~any(isnan(succ_L1_b))
    %                 if weighted
    %                     pU1 = sum(succ_U1_b) / sum(trials_U1_b);
    %                     pL1 = sum(succ_L1_b) / sum(trials_L1_b);
    %                 else
    %                     pU1 = mean(succ_U1_b ./ trials_U1_b);
    %                     pL1 = mean(succ_L1_b ./ trials_L1_b);
    %                 end
    %                 pU1 = max(eps_prop, min(1 - eps_prop, pU1));
    %                 pL1 = max(eps_prop, min(1 - eps_prop, pL1));
    %                 boot_logOR_CT1(i) = log(pU1/(1-pU1)) - log(pL1/(1-pL1));
    %             end
    % 
    %             succ_U2_b = nu2(idx_U2, i);
    %             succ_L2_b = nl2(idx_L2, i);
    %             trials_U2_b = N_U2(idx_U2);
    %             trials_L2_b = N_L2(idx_L2);
    % 
    %             if ~any(isnan(succ_U2_b)) && ~any(isnan(succ_L2_b))
    %                 if weighted
    %                     pU2 = sum(succ_U2_b) / sum(trials_U2_b);
    %                     pL2 = sum(succ_L2_b) / sum(trials_L2_b);
    %                 else
    %                     pU2 = mean(succ_U2_b ./ trials_U2_b);
    %                     pL2 = mean(succ_L2_b ./ trials_L2_b);
    %                 end
    %                 pU2 = max(eps_prop, min(1 - eps_prop, pU2));
    %                 pL2 = max(eps_prop, min(1 - eps_prop, pL2));
    %                 boot_logOR_CT2(i) = log(pU2/(1-pU2)) - log(pL2/(1-pL2));
    %             end
    %         end
    % 
    %         bValid = isfinite(boot_logOR_CT1) & isfinite(boot_logOR_CT2);
    %         if sum(bValid) > 2
    %             boot_r_corr(b) = corr(boot_logOR_CT1(bValid), boot_logOR_CT2(bValid), 'type', 'Pearson');
    %             p_b = polyfit(boot_logOR_CT1(bValid), boot_logOR_CT2(bValid), 1);
    %             y_boot_corr(b, :) = polyval(p_b, x_fit_corr);
    %         else
    %             boot_r_corr(b) = NaN;
    %             y_boot_corr(b, :) = NaN;
    %         end
    %     end
    % 
    %     ci_r_lo = prctile(boot_r_corr, 2.5);
    %     ci_r_hi = prctile(boot_r_corr, 97.5);
    % 
    %     y_reg_lo = prctile(y_boot_corr, 2.5, 1);
    %     y_reg_hi = prctile(y_boot_corr, 97.5, 1);
    %     y_reg_med = median(y_boot_corr, 1, 'omitnan');
    % 
    %     figB = figure('Units', 'centimeters', 'Position', [13 2 16 14], 'Color', 'w');
    %     axB = axes('Parent', figB); hold(axB, 'on');
    % 
    %     fill(axB, [x_fit_corr, fliplr(x_fit_corr)], [y_reg_lo, fliplr(y_reg_hi)], ...
    %         [0.3 0.3 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    %     plot(axB, x_fit_corr, y_reg_med, 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
    % 
    %     for i = 1:numel(xv)
    %         errorbar(axB, xv(i), yv(i), yel(i), yeh(i), xel(i), xeh(i), ...
    %             'color', [0.4 0.4 0.4], 'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 0.9);
    %     end
    %     scatter(axB, xv, yv, 30, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 1);
    % 
    %     for i = 1:numel(xv)
    %         text(axB, xv(i) + 0.05, yv(i) + 0.05, areasPair{i}, 'FontSize', 7, ...
    %             'Color', [0.3 0.3 0.3]);
    %     end
    % 
    %     xlims_B = get(axB, 'XLim'); ylims_B = get(axB, 'YLim');
    %     plot(axB, [0 0], ylims_B, '--k', 'LineWidth', 0.7);
    %     plot(axB, xlims_B, [0 0], '--k', 'LineWidth', 0.7);
    % 
    %     lims = [min([xlims_B(1) ylims_B(1)]), max([xlims_B(2) ylims_B(2)])];
    %     plot(axB, lims, lims, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    % 
    %     xlabel(axB, 'Log-odds (upper/lower) — CT1', 'FontSize', 11);
    %     ylabel(axB, 'Log-odds (upper/lower) — CT2', 'FontSize', 11);
    %     tit = title(axB, sprintf('r = %.2f [%.2f, %.2f], p = %.3g', ...
    %         r_obs, ci_r_lo, ci_r_hi, p_obs));
    %     tit.FontWeight = 'normal'; tit.FontSize = 11;
    %     set(axB, 'FontSize', 10, 'TickDir', 'out', 'Box', 'off');
    %     hold(axB, 'off');
    % 
    %     fprintf('\nLog-odds correlation (CT1 vs CT2):\n');
    %     fprintf('  Observed r: %.3f\n', r_obs);
    %     fprintf('  Bootstrap 95%% CI: [%.3f, %.3f]\n', ci_r_lo, ci_r_hi);
    %     fprintf('  p-value: %.3g\n', p_obs);
    % end

    % ---------------- FIG B: CT1 vs CT2 log-odds correlation ----------------
    
    x_est = Boot_eff_CT1;
    y_est = Boot_eff_CT2;
    
    % MAD of bootstrap distributions
    x_mad = nan(nAreas, 1);
    y_mad = nan(nAreas, 1);
    for i = 1:nAreas
        if ~isempty(Boot_distributions_CT1{i})
            x_mad(i) = mad(Boot_distributions_CT1{i}, 1);
        end
        if ~isempty(Boot_distributions_CT2{i})
            y_mad(i) = mad(Boot_distributions_CT2{i}, 1);
        end
    end
    
    validPair = isfinite(x_est) & isfinite(y_est);
    
    if sum(validPair) > 2
        xv = x_est(validPair);
        yv = y_est(validPair);
        xm = x_mad(validPair);
        ym = y_mad(validPair);
        
        [r_obs, p_obs] = corr(xv, yv, 'type', 'Pearson');
        
        nboot_corr = 10000;
        boot_r_corr = zeros(nboot_corr, 1);
        
        x_fit_corr = linspace(min(xv), max(xv), 100);
        y_boot_corr = zeros(nboot_corr, length(x_fit_corr));
        
        rng(42);
        for b = 1:nboot_corr
            idx_U1 = randi(nU1, nU1, 1);
            idx_L1 = randi(nL1, nL1, 1);
            idx_U2 = randi(nU2, nU2, 1);
            idx_L2 = randi(nL2, nL2, 1);
            
            boot_logOR_CT1 = nan(nAreas, 1);
            boot_logOR_CT2 = nan(nAreas, 1);
            
            for i = 1:nAreas
                succ_U1_b = nu1(idx_U1, i);
                succ_L1_b = nl1(idx_L1, i);
                trials_U1_b = N_U1(idx_U1);
                trials_L1_b = N_L1(idx_L1);
                
                if ~any(isnan(succ_U1_b)) && ~any(isnan(succ_L1_b))
                    if weighted
                        pU1 = sum(succ_U1_b) / sum(trials_U1_b);
                        pL1 = sum(succ_L1_b) / sum(trials_L1_b);
                    else
                        pU1 = mean(succ_U1_b ./ trials_U1_b);
                        pL1 = mean(succ_L1_b ./ trials_L1_b);
                    end
                    pU1 = max(eps_prop, min(1 - eps_prop, pU1));
                    pL1 = max(eps_prop, min(1 - eps_prop, pL1));
                    boot_logOR_CT1(i) = log(pU1/(1-pU1)) - log(pL1/(1-pL1));
                end
                
                succ_U2_b = nu2(idx_U2, i);
                succ_L2_b = nl2(idx_L2, i);
                trials_U2_b = N_U2(idx_U2);
                trials_L2_b = N_L2(idx_L2);
                
                if ~any(isnan(succ_U2_b)) && ~any(isnan(succ_L2_b))
                    if weighted
                        pU2 = sum(succ_U2_b) / sum(trials_U2_b);
                        pL2 = sum(succ_L2_b) / sum(trials_L2_b);
                    else
                        pU2 = mean(succ_U2_b ./ trials_U2_b);
                        pL2 = mean(succ_L2_b ./ trials_L2_b);
                    end
                    pU2 = max(eps_prop, min(1 - eps_prop, pU2));
                    pL2 = max(eps_prop, min(1 - eps_prop, pL2));
                    boot_logOR_CT2(i) = log(pU2/(1-pU2)) - log(pL2/(1-pL2));
                end
            end
            
            bValid = isfinite(boot_logOR_CT1) & isfinite(boot_logOR_CT2);
            if sum(bValid) > 2
                boot_r_corr(b) = corr(boot_logOR_CT1(bValid), boot_logOR_CT2(bValid), 'type', 'Pearson');
                p_b = polyfit(boot_logOR_CT1(bValid), boot_logOR_CT2(bValid), 1);
                y_boot_corr(b, :) = polyval(p_b, x_fit_corr);
            else
                boot_r_corr(b) = NaN;
                y_boot_corr(b, :) = NaN;
            end
        end
        
        ci_r_lo = prctile(boot_r_corr, 2.5);
        ci_r_hi = prctile(boot_r_corr, 97.5);
        
        y_reg_lo = prctile(y_boot_corr, 2.5, 1);
        y_reg_hi = prctile(y_boot_corr, 97.5, 1);
        y_reg_med = median(y_boot_corr, 1, 'omitnan');
        
        figure; set(gcf, 'position', [329.8000 326.6000 482.4000 344.8000]);
        hold on;
        
        fill([x_fit_corr, fliplr(x_fit_corr)], [y_reg_lo, fliplr(y_reg_hi)], ...
            [0.3 0.3 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot(x_fit_corr, y_reg_med, 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
        
        for i = 1:numel(xv)
            errorbar(xv(i), yv(i), ym(i), ym(i), xm(i), xm(i), ...
                'color', [0.4 0.4 0.4], 'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 0.9);
        end
        scatter(xv, yv, 30, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 1);
        
        set(gca, 'FontSize', 12)
        xlabel('log-odds (upper/lower) CT1');
        ylabel('log-odds (upper/lower) CT2');
        
        tit = title(sprintf('r=%.2f [%.2f, %.2f], p=%.3g', r_obs, ci_r_lo, ci_r_hi, p_obs));
        tit.FontWeight = 'normal';
        tit.FontSize = 12;
        
        % Bootstrap distribution
        figure;
        histogram(boot_r_corr, 50, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.6]);
        hold on;
        xline(r_obs, 'r-', 'LineWidth', 2, 'Label', sprintf('r=%.2f', r_obs));
        xline(ci_r_lo, 'b--', 'LineWidth', 1.5);
        xline(ci_r_hi, 'b--', 'LineWidth', 1.5);
        xlabel('Pearson r');
        ylabel('Probability');
        title(sprintf('Bootstrap Distribution: r=%.2f [%.2f, %.2f]', r_obs, ci_r_lo, ci_r_hi));
        set(gca, 'FontSize', 12);
        
        fprintf('\nLog-odds correlation (CT1 vs CT2):\n');
        fprintf('  Observed r: %.3f\n', r_obs);
        fprintf('  Bootstrap mean r: %.3f\n', mean(boot_r_corr, 'omitnan'));
        fprintf('  Bootstrap median r: %.3f\n', median(boot_r_corr, 'omitnan'));
        fprintf('  Bootstrap std: %.3f\n', std(boot_r_corr, 'omitnan'));
    end

    % % ---------------- FIG B: CT1 vs CT2 (zero axes, clear highlights) ----------------
    % % Requirements: Results table with fields
    % %   LME_logOdds_U_vs_L_CT1, LME_CIlo_CT1, LME_CIhi_CT1,
    % %   LME_logOdds_U_vs_L_CT2, LME_CIlo_CT2, LME_CIhi_CT2,
    % %   Area
    % % Optional precomputed visual sizes: sValstot (vector same length as Results)
    % % Colors expected to exist: ColorIntPos, ColorIntNeg, ColorNonSig
    % 
    % % Defensive defaults (if missing)
    % if ~exist('sValstot','var') || isempty(sValstot)
    %     sValstot = repmat(48, height(Results), 1);
    % end
    % if ~exist('ColorIntPos','var') || isempty(ColorIntPos)
    %     ColorIntPos = [0.85 0.35 0.35];
    % end
    % if ~exist('ColorIntNeg','var') || isempty(ColorIntNeg)
    %     ColorIntNeg = [0.18 0.62 0.62];
    % end
    % if ~exist('ColorNonSig','var') || isempty(ColorNonSig)
    %     ColorNonSig = [0.78 0.78 0.78];
    % end
    % 
    % % extract
    % x_est = Results.LME_logOdds_U_vs_L_CT1;
    % y_est = Results.LME_logOdds_U_vs_L_CT2;
    % x_lo  = Results.LME_CIlo_CT1;  x_hi  = Results.LME_CIhi_CT1;
    % y_lo  = Results.LME_CIlo_CT2;  y_hi  = Results.LME_CIhi_CT2;
    % areas = Results.Area;
    % ss_all = sValstot(:);
    % 
    % % Use rows where at least one estimate exists (keeps only meaningful points)
    % maskAny = ~(isnan(x_est) & isnan(y_est));
    % xs = x_est(maskAny); ys = y_est(maskAny);
    % xlo = x_lo(maskAny); xhi = x_hi(maskAny);
    % ylo = y_lo(maskAny); yhi = y_hi(maskAny);
    % areasMask = areas(maskAny);
    % ss = ss_all(maskAny);
    % 
    % n = numel(xs);
    % if n == 0
    %     warning('No data to plot in figB');
    %     return;
    % end
    % 
    % % compute extents (based on all finite limits and points)
    % allvals = [xlo(:); xhi(:); ylo(:); yhi(:); xs(:); ys(:)];
    % xvals=[xlo(:); xhi(:); xs(:);]; yvals= [ylo(:); yhi(:); ys(:)];
    % fvx=xvals(isfinite(xvals)); fvy=yvals(isfinite(yvals));
    % mnx = min(fvx); mxx = max(fvx); mny = min(fvy); mxy = max(fvy);
    % padx = 0.03*(mxx-mnx+eps); pady = 0.03*(mxy-mny+eps);
    % xlims = [mnx-padx, mxx+padx];
    % ylims = [mny-pady, mxy+pady];
    % 
    % % flags: highlight only if BOTH CIs are strictly outside 0
    % ci_out1 = (xlo > 0) | (xhi < 0);   % CT1 CI excludes 0
    % ci_out2 = (ylo > 0) | (yhi < 0);   % CT2 CI excludes 0
    % both_CI_out = ci_out1 & ci_out2;   % highlight criterion
    % 
    % % prepare figure and axis
    % figB = figure('Units','centimeters','Position',[13 2 20 18],'Color','w');
    % axB = axes('Parent',figB); hold(axB,'on');
    % xlim(axB, xlims); ylim(axB, ylims);
    % 
    % % visual params (stronger alphas but still subtle)
    % bgRectAlpha = 0.18;        % non-highlight rectangle fill
    % highlightRectAlpha = 0.70; % highlighted rectangle fill
    % rectFaceColorDefault = ColorNonSig; % very light gray for background rects
    % pointAlpha = 0.7;
    % edgeColorForPoints = [0.12 0.12 0.12]; % outlines for all dots
    % neutralPointFace = ColorNonSig-0.1;
    % 
    % % DRAW ORDER:
    % % --- CI + rectangle logic ---
    % 
    % idxBoth = find(both_CI_out);    % highlight criterion
    % idxRest = find(~both_CI_out);  % everything else
    % 
    % % 1) Non-highlighted: draw CI lines only
    % for k = idxRest(:)'
    %     if all(isfinite([xlo(k), xhi(k), ys(k)]))
    %         plot(axB, [xlo(k) xhi(k)], [ys(k) ys(k)], '-', 'Color', ColorNonSig, 'LineWidth', 1);
    %     end
    %     if all(isfinite([ylo(k), yhi(k), xs(k)]))
    %         plot(axB, [xs(k) xs(k)], [ylo(k) yhi(k)], '-', 'Color', ColorNonSig, 'LineWidth', 1);
    %     end
    % end
    % 
    % % 2) Highlighted: draw rectangles + CI lines
    % for k = idxBoth(:)'
    %     % rectangle coords
    %     lx = xlo(k); hx = xhi(k); ly = ylo(k); hy = yhi(k);
    %     % fallback if missing
    %     epsW = 0.02 * (xlims(2)-xlims(1) + eps);
    %     epsH = 0.02 * (ylims(2)-ylims(1) + eps);
    %     if ~isfinite(lx) || ~isfinite(hx), lx = xs(k)-epsW; hx = xs(k)+epsW; end
    %     if ~isfinite(ly) || ~isfinite(hy), ly = ys(k)-epsH; hy = ys(k)+epsH; end
    %     px = [lx hx hx lx]; py = [ly ly hy hy];
    % 
    %     % color assignment
    %     if xlo(k)>0 && ylo(k)>0
    %         c = ColorNonSig*0.7;
    %     elseif xhi(k)<0 && yhi(k)<0
    %         c = ColorNonSig*0.7;
    %     else
    %         c = ColorNonSig*0.7; % fallback
    %     end
    % 
    %     % rectangle
    %     patch('XData', px, 'YData', py, 'FaceColor', c, 'FaceAlpha', 0.35, ...
    %         'EdgeColor','none','Parent',axB);
    % 
    %     % CI lines on top of rectangle
    %     if all(isfinite([xlo(k), xhi(k), ys(k)]))
    %         plot(axB, [xlo(k) xhi(k)], [ys(k) ys(k)], '-', 'Color', c, 'LineWidth', 1.5);
    %     end
    %     if all(isfinite([ylo(k), yhi(k), xs(k)]))
    %         plot(axB, [xs(k) xs(k)], [ylo(k) yhi(k)], '-', 'Color', c, 'LineWidth', 1.5);
    %     end
    % end
    % 
    % % 3) Scatter dots for all (neutral vs highlighted colored separately)
    % scatter(axB, xs(idxRest), ys(idxRest), ss(idxRest), 'o', ...
    %     'MarkerFaceColor', ColorNonSig*0.7, 'MarkerEdgeColor', edgeColorForPoints, ...
    %     'MarkerFaceAlpha', 0.8, 'LineWidth', 0.7);
    % 
    % for k = idxBoth(:)'
    %     if xlo(k)>0 && ylo(k)>0
    %         c = ColorNonSig*0.7;
    %     elseif xhi(k)<0 && yhi(k)<0
    %         c = ColorNonSig*0.7;
    %     else
    %         c = ColorNonSig*0.7;
    %     end
    %     scatter(axB, xs(k), ys(k), ss(k), 'o', ...
    %         'MarkerFaceColor', c, 'MarkerEdgeColor', edgeColorForPoints, ...
    %         'MarkerFaceAlpha', 0.9, 'LineWidth', 0.9);
    % end
    % % 4) regression line computed only on paired finite xs & ys
    % 
    % % --- Deming regression (errors-in-variables) ---
    % 
    % % Extract finite points
    % mask = isfinite(x_est) & isfinite(y_est) & isfinite(x_lo) & isfinite(x_hi) & isfinite(y_lo) & isfinite(y_hi);
    % x = x_est(mask); y = y_est(mask);
    % 
    % % Approximate per-point SE from CI length (95% CI ≈ estimate ± 1.96*SE)
    % SE_x = (x_hi(mask) - x_lo(mask)) ./ (2*1.96);
    % SE_y = (y_hi(mask) - y_lo(mask)) ./ (2*1.96);
    % 
    % % ---------------- Robust GENERALIZED DEMING with PARAMETRIC bootstrap (no resample) ----------------
    % % Requirements: x,y,SE_x,SE_y (column vectors), axB axis handle
    % % Paste this after plotting points & CI rectangles and before hold(axB,'off').
    % x = x(:); y = y(:); SE_x = SE_x(:); SE_y = SE_y(:);
    % if ~exist('axB','var') || ~isgraphics(axB), error('axB must be an axis handle'); end
    % n = numel(x);
    % if numel(SE_x)~=n || numel(SE_y)~=n, error('SE_x and SE_y must match x,y length'); end
    % 
    % % user-tunable parameters
    % if ~exist('B','var') || isempty(B), B = 3000; end           % bootstrap reps (parametric only)
    % plotLineCI = true;                                            % draw pointwise CI band
    % minSEfactor = 1e-6;                                           % minimal relative SE floor
    % denom_eps = 1e-8;                                             % floor for denom in objective
    % penalty_scale = 1e-8;                                         % penalty for leaving bounds if using fminsearch
    % 
    % % --- diagnostics (print these to help debugging) ---
    % fprintf('DEBUG: n=%d, x span=%.4g, y span=%.4g\n', n, max(x)-min(x), max(y)-min(y));
    % fprintf('DEBUG: SE_x: min=%.4g med=%.4g max=%.4g\n', min(SE_x), median(SE_x), max(SE_x));
    % fprintf('DEBUG: SE_y: min=%.4g med=%.4g max=%.4g\n', min(SE_y), median(SE_y), max(SE_y));
    % 
    % % floor SEs to avoid zeros
    % span_x = max(x)-min(x)+eps; span_y = max(y)-min(y)+eps;
    % minSEx = max(min(SE_x(SE_x>0))/10, minSEfactor*span_x);
    % minSEy = max(min(SE_y(SE_y>0))/10, minSEfactor*span_y);
    % SE_x = max(SE_x, minSEx);
    % SE_y = max(SE_y, minSEy);
    % 
    % % --- choose sensible slope bounds based on data ---
    % % naive slope from ranges (avoid divide by zero)
    % % --- robust slope-guess & bound computation (replace previous Sbound code) ---
    % % safe slope from OLS (scalar)
    % if var(x) > 0
    %     % robust scalar OLS slope: cov(x,y)/var(x) computed in scalar form
    %     slope_ols = sum( (x - mean(x)) .* (y - mean(y)) ) / ( sum( (x - mean(x)).^2 ) + eps );
    % else
    %     slope_ols = 0;
    % end
    % 
    % % slope from range heuristic (avoid divide-by-zero)
    % slope_range_guess = ( max(y) - min(y) ) / ( max(x) - min(x) + eps );
    % 
    % % choose slope_dem_guess: prefer slope_ols (scalar)
    % slope_dem_guess = slope_ols;
    % if ~isfinite(slope_dem_guess)
    %     slope_dem_guess = 0;
    % end
    % 
    % % Sbound: compute from a vector using max(...) on a vector, then cap it
    % Sbound = max([50, 10*abs(slope_dem_guess), 10*abs(slope_range_guess)]);
    % Sbound = min(Sbound, 1e3);  % hard ceiling to avoid ridiculous bounds
    % 
    % fprintf('DEBUG: slope_ols=%.4g, slope_range=%.4g, Sbound=%.4g\n', slope_ols, slope_range_guess, Sbound);
    % 
    % 
    % b1_min = -Sbound; b1_max = Sbound;
    % 
    % % --- stable objective with denominator floor ---
    % objfun = @(b,xv,yv,sex,sey) sum( ((yv - (b(1) + b(2).*xv)).^2) ./ max(sey.^2 + (b(2).^2).*sex.^2, denom_eps) );
    % 
    % % If fmincon available, use it with bounds; otherwise use fminsearch with penalty
    % useFmincon = exist('fmincon','file')==2;
    % 
    % opts_fs = optimset('Display','off','TolX',1e-10,'TolFun',1e-10,'MaxIter',5e3,'MaxFunEvals',5e3);
    % 
    % % helper that enforces bounds when fminsearch is used (soft penalty)
    % if ~useFmincon
    %     objfun_bounded = @(b,xv,yv,sex,sey) objfun(b,xv,yv,sex,sey) + penalty_scale * ( max(0,b(2)-b1_max)^2 + max(0,b1_min-b(2))^2 );
    % end
    % 
    % % initialization for fit on observed data: try multiple inits
    % if var(x) > 0
    %     b1_init = slope_ols;
    % else
    %     b1_init = 0;
    % end
    % b0_init = mean(y) - b1_init*mean(x);
    % inits = {[b0_init; b1_init], [mean(y); 0], [median(y); 0.5*b1_init]};
    % 
    % bestb = []; bestval = inf;
    % for t = 1:numel(inits)
    %     try
    %         if useFmincon
    %             % bounds on [b0, b1] only on b1
    %             A = []; b = []; Aeq = []; beq = []; lb = [-inf; b1_min]; ub = [inf; b1_max];
    %             btry = fmincon(@(bb) objfun(bb,x,y,SE_x,SE_y), inits{t}, A,b,Aeq,beq,lb,ub,[], optimoptions('fmincon','Display','off','Algorithm','interior-point'));
    %             val = objfun(btry,x,y,SE_x,SE_y);
    %         else
    %             btry = fminsearch(@(bb) objfun_bounded(bb,x,y,SE_x,SE_y), inits{t}, opts_fs);
    %             % clip slope to bounds if slight drift occurred
    %             btry(2) = min(max(btry(2), b1_min), b1_max);
    %             val = objfun(btry,x,y,SE_x,SE_y);
    %         end
    %         if isfinite(val) && val < bestval
    %             bestval = val; bestb = btry;
    %         end
    %     catch ME
    %         % ignore, try next init
    %     end
    % end
    % 
    % if isempty(bestb)
    %     error('Generalized Deming fit failed on observed data with bounded optimization. Check SEs.');
    % end
    % beta0_gdem = bestb(1); beta1_gdem = bestb(2);
    % fprintf('DEBUG: fitted gDem slope = %.6g\n', beta1_gdem);
    % 
    % % plot fitted line
    % xgrid = linspace(min(x)-eps, max(x)+eps, 200);
    % y_gdem = beta0_gdem + beta1_gdem * xgrid;
    % hold(axB,'on');
    % h_gdem = plot(axB, xgrid, y_gdem, '-', 'LineWidth', 2, 'Color', [0 0.55 0.2]);
    % 
    % % ---------------- PARAMETRIC BOOTSTRAP (no resampling) ----------------
    % % For each rep: simulate xss* ~ N(x,SE_x^2), yss* ~ N(y,SE_y^2), fit bounded generalized Deming
    % bs_b1 = nan(B,1);
    % bs_b0 = nan(B,1);
    % bs_r  = nan(B,1);
    % failed = 0;
    % 
    % rng(2); % reproducible
    % for ib = 1:B
    %     xss = x + randn(n,1).*SE_x;
    %     yss = y + randn(n,1).*SE_y;
    % 
    %     % fit with same bounded approach
    %     found = false;
    %     for t = 1:numel(inits)
    %         try
    %             if useFmincon
    %                 b_try = fmincon(@(bb) objfun(bb,xss,yss,SE_x,SE_y), inits{t}, [],[],[],[], [-inf;b1_min], [inf;b1_max], [], optimoptions('fmincon','Display','off','Algorithm','interior-point'));
    %             else
    %                 b_try = fminsearch(@(bb) objfun_bounded(bb,xss,yss,SE_x,SE_y), inits{t}, opts_fs);
    %                 b_try(2) = min(max(b_try(2), b1_min), b1_max);
    %             end
    %             % validate b_try
    %             if isfinite(b_try(2)) && abs(b_try(2)) <= 1e6
    %                 bs_b0(ib) = b_try(1);
    %                 bs_b1(ib) = b_try(2);
    %                 found = true;
    %                 break;
    %             end
    %         catch
    %             % try next init
    %         end
    %     end
    %     if ~found
    %         % fallback: OLS on simulated if bounded fit failed
    %         if var(xss) > 0
    %             b1_ols = (cov(xss,yss)/var(xss)); if numel(b1_ols)>1, b1_ols = b1_ols(1,2)/var(xss); end
    %             bs_b1(ib) = min(max(b1_ols, b1_min), b1_max);
    %             bs_b0(ib) = mean(yss) - bs_b1(ib)*mean(xss);
    %         else
    %             bs_b1(ib) = NaN; bs_b0(ib)=NaN;
    %         end
    %         failed = failed + 1;
    %     end
    % 
    %     % correlation
    %     if var(xss) > 0 && var(yss) > 0
    %         rtmp = corr(xss, yss, 'Rows','complete', 'Type','Pearson');
    %         if isfinite(rtmp)
    %             bs_r(ib) = max(min(rtmp, 0.999999), -0.999999);
    %         else
    %             bs_r(ib) = NaN;
    %         end
    %     else
    %         bs_r(ib) = NaN;
    %     end
    % end
    % 
    % % filter finite and reasonable slopes
    % validSlopeMask = isfinite(bs_b1) & abs(bs_b1) < 1e6;
    % finiteSlopes = bs_b1(validSlopeMask);
    % finiteR = bs_r(isfinite(bs_r));
    % 
    % % diagnostics
    % fprintf('DEBUG: bootstrap reps B=%d; failed fits (fell back or failed) = %d (%.1f%%)\n', B, failed, 100*failed/B);
    % if numel(finiteSlopes) < max(10, 0.1*B)
    %     warning('Too few successful bootstrap slope fits (%d). Results may still be unstable.', numel(finiteSlopes));
    % end
    % 
    % % compute slope CI (percentile)
    % if isempty(finiteSlopes)
    %     slope_median = beta1_gdem;
    %     slope_CI = [NaN NaN];
    % else
    %     slope_median = median(finiteSlopes);
    %     slope_CI = prctile(finiteSlopes, [2.5 97.5]);
    % end
    % 
    % % % compute r CI
    % % if isempty(finiteR)
    % %     r_med = NaN; r_CI = [NaN NaN];
    % % else
    % %     r_med = median(finiteR);
    % %     r_CI = prctile(finiteR, [2.5 97.5]);
    % % end
    % 
    % % optional: plot pointwise CI band from bootstrap fits
    % if plotLineCI && ~isempty(finiteSlopes)
    %     validIdx = validSlopeMask;
    %     preds = bs_b0(validIdx) + bs_b1(validIdx) .* xgrid;
    %     y_lo = prctile(preds, 2.5, 1);
    %     y_hi = prctile(preds, 97.5, 1);
    %     patch(axB, [xgrid fliplr(xgrid)], [y_lo fliplr(y_hi)], [0 0.55 0.2], 'FaceAlpha', 0.16, 'EdgeColor', 'none');
    %     uistack(h_gdem,'top');
    % end
    % 
    % % 3. Pearson r with Fisher-z CI
    % [r,p] = corr(x,y,'Type','Pearson','Rows','complete');
    % z = atanh(r);
    % n = numel(x);
    % se = 1/sqrt(n-3);
    % zCI = z + se*[-1 1]*1.96;
    % rCI = tanh(zCI);
    % 
    % % annotate axis with robust stats
    % slopeStr = sprintf('gDem slope = %.3g [%.3g, %.3g]', beta1_gdem, slope_CI(1), slope_CI(2));
    % % rStr = sprintf('r = %.2f [%.2f, %.2f]', r_med, r_CI(1), r_CI(2));
    % rStr=sprintf('Pearson r = %.2f [%.2f, %.2f]', r, rCI(1), rCI(2));
    % text(axB, 0.02, 0.98, sprintf('%s\n%s', slopeStr, rStr), 'Units','normalized','VerticalAlignment','top', ...
    %     'BackgroundColor',[1 1 1 0.9], 'EdgeColor',[0.85 0.85 0.85], 'FontSize',10);
    % 
    % % store bootstrap arrays in base for inspection (optional)
    % assignin('base','gDem_bs_b1', bs_b1);
    % assignin('base','gDem_bs_r', bs_r);
    % assignin('base','gDem_bs_b0', bs_b0);

    % 
    % % hold(axB,'off');
    % %
    % % paired = isfinite(xs) & isfinite(ys);
    % % if sum(paired) >= 2
    % %     x_p = xs(paired); y_p = ys(paired);
    % %     % linear fit
    % %     b = polyfit(x_p, y_p, 1);
    % %     xl = linspace(min(x_p)-eps, max(x_p)+eps, 200);
    % %     yl = polyval(b, xl);
    % %     % draw regression line (subtle color)
    % %     regColor = [0.45 0.35 0.35]; % neutral muted
    % %     hreg = plot(axB, xl, yl, '-', 'LineWidth', 1.6, 'Color', regColor);
    % %     % compute r and p and slope
    % %     [r, p_r] = corr(x_p, y_p, 'Rows','pairwise', 'Type','Pearson');
    % %     slope = b(1);
    % %     statsStr = sprintf('r = %.2f (p = %.4g), slope = %.2f', r, p_r, slope);
    % % else
    % %     statsStr = 'r = NaN';
    % % end
    % % 
    % % 
    % 
    % % 6) labels for highlighted items (deterministic placement: top-right inside rectangle when possible)
    % text_dx = 0.01 * (xlims(2)-xlims(1));
    % text_dy = 0.01 * (ylims(2)-ylims(1));
    % for ii = 1:numel(idxBoth)
    %     k = idxBoth(ii);
    %     % place labe
    %     tx = xs(k); ty = ys(k); 
    %     if xs(k)>0 && ys(k)>0 
    %         hal = 'left'; val = 'bottom';
    %     elseif xs(k)>0 && ys(k)<0 
    %         hal = 'left'; val = 'bottom';
    %     elseif xs(k)<0 && ys(k)<0
    %         hal = 'right'; val = 'top';
    %     else xs(k)<0 && ys(k)>0;
    %         hal = 'right'; val = 'top';
    %     end
    %     txt = char(areasMask(k));
    %     text(axB, tx, ty, txt, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', hal, ...
    %         'VerticalAlignment', val, 'Interpreter', 'none', 'Color', 'k', 'BackgroundColor','none');
    % end
    % 
    % % 7) 0 axes dotted on top
    % xlims_now = get(axB,'XLim'); ylims_now = get(axB,'YLim');
    % plot(axB, [0 0], ylims_now, '--', 'Color', [0 0 0], 'LineWidth', 0.7);
    % plot(axB, xlims_now, [0 0], '--', 'Color', [0 0 0], 'LineWidth', 0.7);
    % 
    % % % 8) print stats top-left normalized
    % % text(axB, 0.02, 0.98, statsStr, 'Units','normalized', 'VerticalAlignment','top', 'FontSize',9);
    % 
    % % tidy axis labels and title
    % xlabel(axB,'LME log-odds (U vs L) — CT1');
    % ylabel(axB,'LME log-odds (U vs L) — CT2');
    % title(axB,'CT1 vs CT2: effects per area (highlight = both CI outside 0)', 'FontWeight','bold');
    % set(axB, 'TickDir','out', 'FontSize', 10, 'Box', 'off');
    % hold(axB,'off');

end

end  % end of main function

% -------------------- Local helper validators --------------------
function ok = validateColor(x)
ok = isnumeric(x) && numel(x)==3 && all(x>=0) && all(x<=1);
if ~ok, error('Color must be 1x3 RGB (0..1).'); end
end

function ok = validatePalette(x)
ok = isnumeric(x) && size(x,2)==3 && all(x(:)>=0) && all(x(:)<=1);
if ~ok, error('Palette must be Mx3 RGB (0..1).'); end
end

function draw_forest_aligned(ax, effects, ci_lo, ci_hi, names, qvals, colPos, colNeg, nonSigCol, svals, titleStr, orderRows, topN)
if nargin < 13, topN = 0; end
axes(ax); hold(ax,'on');

nTotal = numel(names);
if nargin < 11 || isempty(orderRows)
    orderRows = 1:nTotal;
end

ypos = 1:nTotal;

allVals = [ci_lo(:); ci_hi(:); effects(:)];
finiteVals = allVals(isfinite(allVals));
if isempty(finiteVals)
    minX = -1; maxX = 1;
else
    minX = min(finiteVals); maxX = max(finiteVals);
    rngSpan = (maxX - minX) + eps;
    pad = 0.08*rngSpan;
    minX = minX - pad; maxX = maxX + pad;
    maxAbs = max(abs([minX, maxX]));
    minX = -1.05*maxAbs; maxX = 1.05*maxAbs;
end

for r = 1:nTotal
    idx = orderRows(r);
    e = effects(idx);
    l = ci_lo(idx);
    h = ci_hi(idx);
    qk = qvals(idx);
    sv = svals(idx);

    if ~isfinite(e) && ~isfinite(l) && ~isfinite(h)
        continue;
    end

    % if ~isnan(qk) && (qk < 0.05)
    %     if e > 0, col = colPos; else col = colNeg; end
    % else
    %     col = nonSigCol;
    % end

    if e > 0.4 && l > 0
        col = colPos; 
    elseif e < -0.4 && h < 0
        col = colNeg; 
    else
        col = nonSigCol;
    end
    lineCol = max(0, col);

    if isfinite(l) && isfinite(h)
        plot(ax, [l h], [ypos(r) ypos(r)], '-', 'LineWidth', 1.6, 'Color', [lineCol 0.6]);
    end
    if isfinite(e)
        scatter(ax, e, ypos(r), sv, 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', lineCol, 'MarkerFaceAlpha',0.6);
    end
end

plot(ax, [0 0], [0 nTotal+1], '--k', 'LineWidth', 1);

yticklabels_ordered = names(orderRows);
set(ax, 'YTick', ypos, 'YTickLabel', yticklabels_ordered, 'YDir', 'reverse', 'FontSize', 9, 'TickDir', 'out');

xlim(ax, [minX maxX]);
xlabel(ax, 'log-odds (upper vs lower)');
title(ax, titleStr, 'FontWeight', 'bold');

box(ax,'off');
hold(ax,'off');
end
