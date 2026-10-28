function N = group_totals(totalIn, counts, nAnim, name)
%GROUP_TOTALS Total inputs per animal as a column vector.
%   N = GROUP_TOTALS(totalIn, counts, nAnim, name) returns totalIn as a
%   column (a scalar is repeated for all nAnim animals). If totalIn is
%   empty, the row sums of counts are used and a warning names the input.

if isempty(totalIn)
    N = sum(counts, 2);
    warning('%s not provided: defaulting to row sums of the counts.', name);
elseif isscalar(totalIn)
    N = repmat(totalIn, nAnim, 1);
else
    N = totalIn(:);
end
assert(numel(N) == nAnim, '%s length mismatch.', name);
end
