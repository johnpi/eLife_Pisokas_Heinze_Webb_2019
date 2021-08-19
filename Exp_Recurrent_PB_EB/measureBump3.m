%% Another version of the bump measurement
% Detect if there was a bump at the end of the simulation
% Measure activity bump width as width at 50% of maximum
function bump=measureBump3(Activity, fixedThreshold, weightVectorsOnProximityToMax)
    N = size(Activity, 2);    % Number of neurons
    neuron_idxs = [1:N];      % The indeces of the neurons we are considering
    % Get the angular azimuth that each neuron corresponds to
    theta_angles = getOrientation(neuron_idxs, N);
    
    vectors = add_vectors(Activity, theta_angles, weightVectorsOnProximityToMax);
    
    x_total_list     = vectors(:,1);
    y_total_list     = vectors(:,2);
    r_total_list     = vectors(:,3);
    theta_total_list = vectors(:,4);
    
    bump.r_ts = r_total_list;
    bump.theta_ts = theta_total_list;
end

function orientations=getOrientation(idx_list, N)
    % idx_list : list indeces of neurons
    % N        : Number of neurons
    % returns  : list of theta angle preferences of neurons in the list of indeces
    % 
    % The circle 360 is broken into N intervals and the angle corresponding in 
    % the middle of the interval is calculated. 
    interval = 360 / N;
    interval_middle = 360 / N / 2;
    
    orientations = interval_middle + interval * (idx_list - 1); % Make the idx_list 0 based
end

function vectors=add_vectors(rs_list, theta_angles, weightVectorsOnProximityToMax)
    % rs_list an array with each row containing the activity values of N
    % neurons at one time step. Each row corresponds to one time step. 
    % Each columns corresponds to one neuron. 
    % theta_angles is a row vector containing the theta agnle corresponding
    % to each column (neuron) of rs_list. 
    theta_rads = deg2rad(theta_angles);
    % Make sure larger than 360deg angles wrap around
    theta_rads = mod(theta_rads, 2*pi);
    
    if ~weightVectorsOnProximityToMax
        % Convert the vectors from polar coordinates to cartesian coordinates
        x = rs_list .* cos(theta_rads);
        y = rs_list .* sin(theta_rads);
    else
        % This option attempts to filter out the effect of spurious noise
        % to the calculation of the heading vector by giving more
        % importance to the longest vector and its neighbours. This avoids
        % the effect of instantaneous noise say 180deg away from the bump 
        % peak moving the average heading in the middle between the actual
        % heading and the noise. 
        % Calculate the weighted vectors
        % Find the position of the max element in each row (longest vector)
        [max_val, max_idx] = max(rs_list, [], 2);

        % Calculate the weighting of each vector in each row with the vector 
        % with maximum length having weight 1 and the other less proportional 
        % to their angular distance from the longest vector following a 
        % von Mises distribution. 
        kappa = pi*2; % Width of von Mises distribution. kappa=pi*2 results to Gaussian with 61deg FWHM
        kappa = pi/100; % Width of von Mises distribution. kappa=pi*2 results to Gaussian with 61deg FWHM
        baseline = 0.5;
        [weights angles] = circ_vmpdf(theta_rads, theta_rads(max_idx), kappa);
        % Normalise p to the range of <baseline> to 1 because it varies with kappa
        weights = (weights - min(min(weights))) * (1 - baseline) / max(max(weights - min(min(weights)))) + baseline;

        % Convert the weighted vectors from polar coordinates to cartesian coordinates
        x = rs_list .* weights' .* cos(theta_rads);
        y = rs_list .* weights' .* sin(theta_rads);
    end
    
    % Calculate the total x and y coordinates of the the vectors described
    % in each row (time step).
    x_total_list = sum(x, 2);
    y_total_list = sum(y, 2);
    
    % Convert cartesian coordinates to polar coordinates
    r_total_list = sqrt(x_total_list.^2 + y_total_list.^2);
    theta_total_list = atan2(y_total_list, x_total_list); % This returns (-pi,pi]
    % Convert the negative angle values to positive so the range is [0,2*pi)
    neg_elements = y_total_list < 0;
    theta_total_list(neg_elements) = theta_total_list(neg_elements) + 2*pi;
    
    % Convert to degrees
    theta_total_list = rad2deg(theta_total_list);
    
    vectors = [x_total_list y_total_list r_total_list theta_total_list];
end
