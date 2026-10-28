function prop = group_proportion(succ, trials, weighted)
%GROUP_PROPORTION Input fraction of a group of animals.
%   prop = GROUP_PROPORTION(succ, trials, weighted) returns the mean of the
%   per-animal fractions succ./trials (weighted = false) or the pooled
%   fraction sum(succ)/sum(trials) (weighted = true).

if weighted
    prop = sum(succ) / sum(trials);
else
    prop = mean(succ ./ trials);
end
end
