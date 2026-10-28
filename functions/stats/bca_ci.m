function [lo, hi] = bca_ci(bootDist, z0, a, alphaLevel)
%BCA_CI Bias-corrected and accelerated bootstrap interval.
%   [lo, hi] = BCA_CI(bootDist, z0, a, alphaLevel) returns the
%   (1 - alphaLevel) BCa interval of the bootstrap distribution bootDist.
%   The adjusted percentiles are limited to [0.1, 99.9].

zl = norminv(alphaLevel / 2);
zu = norminv(1 - alphaLevel / 2);
al = normcdf(z0 + (z0 + zl) / (1 - a * (z0 + zl)));
au = normcdf(z0 + (z0 + zu) / (1 - a * (z0 + zu)));
al = max(0.001, min(0.999, al));
au = max(0.001, min(0.999, au));
ci = prctile(bootDist, [al * 100, au * 100]);
lo = ci(1);
hi = ci(2);
end
