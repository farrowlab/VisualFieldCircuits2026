function z0 = bca_bias_correction(bootDist, theta_hat, nBootstrap)
%BCA_BIAS_CORRECTION Bias-correction constant z0 of the BCa interval.
%   z0 = norminv(fraction of bootstrap estimates below theta_hat), limited
%   to [-3, 3] when the fraction is 0 or 1.

prop_less = sum(bootDist < theta_hat) / nBootstrap;
if prop_less == 0
    z0 = -3;
elseif prop_less == 1
    z0 = 3;
else
    z0 = norminv(prop_less);
end
end
