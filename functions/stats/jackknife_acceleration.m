function a = jackknife_acceleration(jack)
%JACKKNIFE_ACCELERATION Acceleration constant of the BCa interval.
%   a = JACKKNIFE_ACCELERATION(jack) from the leave-one-out estimates jack;
%   0 if undefined or |a| > 1.

jack_mean = mean(jack);
numerator = sum((jack_mean - jack).^3);
denominator = 6 * (sum((jack_mean - jack).^2))^(3/2);
if denominator == 0 || isnan(denominator) || abs(numerator / denominator) > 1
    a = 0;
else
    a = numerator / denominator;
end
end
