function [y_normpdf_normalised]=weights_norm_dist(sigma, num_of_samples)
    sample_points = pi / (num_of_samples-1);
    t = 0:sample_points:pi;
    y_norm=normpdf(t, pi/2, sigma);
    y_normpdf_normalised = y_norm/max(y_norm);
end