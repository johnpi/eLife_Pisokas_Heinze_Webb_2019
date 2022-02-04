%% 
%% inhibition_distr_type : one of tom_sin_approx or gaussian
%% num_of_samples        : length of the returned vector
%% inhibition_width_sigma: if using the Gaussian distribution, the sigma
%% 
function [inhibition_dist]=get_inhib_distr_vector(inhibition_distr_type, num_of_samples, inhibition_width_sigma)
    if strcmp(inhibition_distr_type, 'uniform')
        inhibition_dist = ones(1, num_of_samples);
    elseif strcmp(inhibition_distr_type, 'tom_sin_approx')
        inhibition_dist = [0 0.14644661 0.5 0.85355339 1.0 0.85355339 0.5 0.14644661 0];
    elseif strcmp(inhibition_distr_type, 'gaussian')
        % Get values for the inhibition strengths using gaussian distribution
        if exist('inhibition_width_sigma', 'var')
            sigma = inhibition_width_sigma;
        else
            sigma = 0.5;        
        end
        inhibition_dist = weights_norm_dist(sigma, num_of_samples);
    end

end