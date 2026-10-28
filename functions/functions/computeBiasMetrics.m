function M = computeBiasMetrics(pU, pL, varargin)
% computeBiasMetrics  Compute several bias metrics between two fractions.
%   M = computeBiasMetrics(pU, pL) where pU/pL are vectors (same length),
%   either in 0..1 or 0..100 (function auto-detects).
%
% Optional name/value:
%   'Alpha'   : pseudocount added to each fraction (default: auto)
%   'TotalN'  : effective total input neurons (scalar) -- used to set Alpha
%   'Plot'    : true/false (default false) - make a diagnostic plot
%
% Returns table M with columns:
%   pU, pL, raw_log2, reg_log2, mean_frac, wb_mean, wb_sum, wb_sqrt, logit_diff_w

% parse inputs
p = inputParser;
addRequired(p,'pU',@(x)isnumeric(x));
addRequired(p,'pL',@(x)isnumeric(x));
addParameter(p,'Alpha',[],@(x)isnumeric(x) && isscalar(x));
addParameter(p,'TotalN',[],@(x)isnumeric(x) && isscalar(x));
addParameter(p,'Plot',false,@islogical);
parse(p,pU,pL,varargin{:});

pU = p.Results.pU(:);
pL = p.Results.pL(:);
assert(numel(pU)==numel(pL),'pU and pL must be same length');

% detect scale: if values > 1 assume percentages 0..100
if max([pU; pL]) > 1
    pU = pU ./ 100;
    pL = pL ./ 100;
end

% small epsilon to avoid log(0)
epsv = 1e-12;

% choose alpha / pseudocount
if ~isempty(p.Results.Alpha)
    alpha = p.Results.Alpha;
elseif ~isempty(p.Results.TotalN)
    % if you have an effective total N, use pseudocount = 1/(2*N) (Jeffreys-ish)
    alpha = 1 / (2 * p.Results.TotalN);
else
    % fallback: set alpha to half of the smallest non-zero fraction (or a tiny floor)
    nonzero = min([pU(pU>0); pL(pL>0)]);
    if isempty(nonzero)
        alpha = 1e-6;
    else
        alpha = max(nonzero/2, 1e-8);
    end
end

% raw metrics
raw_log2 = log2( (pU + epsv) ./ (pL + epsv) );    % plain log2 ratio
reg_log2 = log2( (pU + alpha) ./ (pL + alpha) );  % regularized

mean_frac = (pU + pL)/2;
sum_frac = pU + pL;

wb_mean = mean_frac .* raw_log2;   % recommend this
wb_sum  = sum_frac  .* raw_log2;
wb_sqrt = sqrt(mean_frac) .* raw_log2;

% logit difference (use stabilized logit)
logit = @(x) log((x + epsv)./(1 - x + epsv));
logit_diff = logit(pU) - logit(pL);
logit_diff_w = mean_frac .* logit_diff;

% put into table
M = table(pU, pL, mean_frac, raw_log2, reg_log2, wb_mean, wb_sum, wb_sqrt, logit_diff, logit_diff_w);

% optional plot

figure('Units','centimeters','Position',[2 2 18 10]);
subplot(1,2,1);
scatter(raw_log2, mean_frac, 36, 'filled');
xlabel('log_2(p_U / p_L)'); ylabel('mean fraction'); title('raw log2 vs abundance');
subplot(1,2,2);
scatter(wb_mean, mean_frac, 36, 'filled');
xlabel('weighted bias (mean * log2)'); ylabel('mean fraction'); title('weighted vs abundance');


end