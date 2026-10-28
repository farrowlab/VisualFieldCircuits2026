function q = bh_qvalues(pvals)
%BH_QVALUES Benjamini-Hochberg q-values (column); NaN p-values stay NaN.

q = nan(numel(pvals), 1);
valid = ~isnan(pvals);
if any(valid)
    q(valid) = mafdr(pvals(valid), 'BHFDR', true);
end
end
