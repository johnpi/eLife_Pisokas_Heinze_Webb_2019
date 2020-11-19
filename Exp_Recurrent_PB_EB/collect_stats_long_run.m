%% This script runs and plots multiple simulations 
%% running the simulation for long time. It collects 
%% the spiking activity over time to see how the bump 
%% drifts and deviates. It also takes as parameters and 
%% changes specific synaptic type by specified percentage 
%% to study the effect of noise in the circuit.

%% ARGUMENTS:
%% exp_model : Set it to 'dros' or 'locust' to run a drosophila or 
%%             locust simulation. Other available options 'mixed_species' 
%%             'dros_without_P_EG' and 'locust_without_P_EG'
%% changed_synapse_type : String specifying which synapse type to add noise 
%%                 to. Options are: 
%%                                 'none' 
%%                                 'All'
%%                                 'P_EN_to_E_PG_all' 
%%                                 'P_EN_to_E_PG_alleq' 
%%                                 'P_EN_to_E_PG_left' 
%%                                 'P_EN_to_E_PG_right' 
%%                                 'P_EG_to_E_PG_all' 
%%                                 'P_EG_to_E_PG_alleq' 
%%                                 'P_EG_to_E_PG_left' 
%%                                 'P_EG_to_E_PG_right' 
%%                                 'E_PG_to_P_EN_all' 
%%                                 'E_PG_to_P_EN_alleq' 
%%                                 'E_PG_to_P_EN_left' 
%%                                 'E_PG_to_P_EN_right' 
%%                                 'E_PG_to_P_EG_all' 
%%                                 'E_PG_to_P_EG_alleq' 
%%                                 'E_PG_to_P_EG_left' 
%%                                 'E_PG_to_P_EG_right' 
%%                                 'E_PG_to_Pintr_all' 
%%                                 'E_PG_to_Pintr_alleq' 
%%                                 'E_PG_to_Pintr_left' 
%%                                 'E_PG_to_Pintr_right' 
%%                                 'Pintr_to_P_EN_all' 
%%                                 'Pintr_to_P_EG_all' 
%%                                 'Pintr_to_Pintr_all' 
%%                                 'Pintr_to_P_EN_alleq' 
%%                                 'Pintr_to_P_EG_alleq' 
%%                                 'Pintr_to_Pintr_alleq' 
%% percent_noise_values : Vector with values of noise percentage to add 
%%                      to the synaptic weights. Eg -50 or [-50]
%%                      If this argument is equal to [-1, -1] its values are 
%%                      ignored and the default noise values are used.
%% display_plots : Set 0 or 1 to not or do display a plot of the spiking 
%%                 activity at the end of each run. 
%% data_filename : A filename where to store the collected data. Eg
%%                 'Drift_Analysis/Data/collect_stats_long_run_data.mat'
%% stim_type     : 'init_only' ?init_and_jump_bump? 
%%                 'init_and_unihemispheric' ?init_and_unihemispheric?
%%                 ?init_only_vonMises? or ?real_trajectoy_vonMises?
%% stim_neuron_list : vector of neuron numbers to apply the initial 
%%                 stimulus to. Numbers are in reference to the inputList 
%%                 range of input neurons. So a stim_neuron_list = [1,2]
%%                 will result to stimulus being applied to the first and 
%%                 second neurons in inputList.
%%
%% OUTPUT FILE FORMAT:
%% The output is writen in a file as an array of records
%%     [RECORD1; RECORD2; ... RECORDn]
%% each record is a cell array with these fields
%%        data_record{1}  = exp_model; ('dros' or 'locust')
%%        data_record{2}  = sparse(logical(spikeSequence));
%%        data_record{3}  = n_P_EN;
%%        data_record{4}  = n_P_EG;
%%        data_record{5}  = n_E_PG;
%%        data_record{6}  = n_Pintr;
%%        data_record{7}  = con_matrix_filename;
%%        data_record{8}  = con_matrix_parameter_set; ([... ... ...])
%%        data_record{9}  = inhibition_distr_type;
%%        data_record{10} = inhibition_width_sigma;
%%        data_record{11} = dt; (time step of simulation)
%%        data_record{12} = Simulation Duration;
%%        data_record{13} = changed_synapse_type;
%%        data_record{14} = change_percent;
%%        data_record{15} = con_matrix;
%%        Further data_record entries are added depending on the length of 
%%        additional stimulus parameters.
%%

%% Example runs:
%%
%%    Two neighbouring neurons stimulated (4 and 5) (neurons numbered 1-8 or 1-9)
%%    for i=[1:10] disp(sprintf('Iteration %d\n', i)); j=4; stim_j = [j, mod((j+1)-1, 9)+1]; stim_strength_fact=1; collect_stats_long_run('dros',   'none',              [0],      1, 'Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_dros_2_0.mat', 'init_only', stim_j, stim_strength_fact) ; end
%%    for i=[1:10] disp(sprintf('Iteration %d\n', i)); j=4; stim_j = [j, mod((j+1)-1, 8)+1]; stim_strength_fact=1; collect_stats_long_run('locust', 'none',              [0],      1, 'Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_locust_2_0.mat', 'init_only', stim_j, stim_strength_fact) ; end

%%    Stimulating group of neurons around a single centre with von Mises distribution.
%%    for j=[1:10] for i=[1:0.25:9] disp(sprintf('Iteration stimulate E-PG neuron %d with von Mises\n', i)); stim_j = [i]; stim_strength_fact=1; collect_stats_long_run('dros',   'none',  [0], 1, 'Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_dros_6_0.mat', 'init_only_vonMises', stim_j, stim_strength_fact) ; end ; end
%%    for j=[1:10] for i=[1:0.25:8] disp(sprintf('Iteration stimulate E-PG neuron %d with von Mises\n', i)); stim_j = [i]; stim_strength_fact=1; collect_stats_long_run('locust', 'none',  [0], 1, 'Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_locust_6_0.mat', 'init_only_vonMises', stim_j, stim_strength_fact) ; end ; end

%% Bump jumping with the spiking model or with the rate based model
%%    for i = [1:10] disp(sprintf('Iteration = %d\n', i)); j=2; stim_j = [j]; stim_strength_fact=1; collect_stats_long_run('dros',     'none',  [0], 1, 'Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_dros_jump_spiking_0.mat',   'init_and_jump_bump', stim_j, stim_strength_fact, 4.0, 'spiking', 180); end
%%    for i = [1:10] disp(sprintf('Iteration = %d\n', i)); j=2; stim_j = [j]; stim_strength_fact=1; collect_stats_long_run('dros',     'none',  [0], 1, 'Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_dros_jump_ratebased_0.mat',   'init_and_jump_bump', stim_j, stim_strength_fact, 4.0, 'ratebased', 180); end
%%    for i = [1:10] disp(sprintf('Iteration = %d\n', i)); j=2; stim_j = [j]; stim_strength_fact=1; collect_stats_long_run('locust',   'none',  [0], 1, 'Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_locust_jump_spiking_0.mat', 'init_and_jump_bump', stim_j, stim_strength_fact, 4.0, 'spiking', 180); end
%%    for i = [1:10] disp(sprintf('Iteration = %d\n', i)); j=2; stim_j = [j]; stim_strength_fact=1; collect_stats_long_run('locust',   'none',  [0], 1, 'Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_locust_jump_ratebased_0.mat', 'init_and_jump_bump', stim_j, stim_strength_fact, 4.0, 'ratebased', 180); end

%% Bump jumping with random noise in the synaptic weights of all synapse types
%%    trials = 10; noise_on_synapse_type = 'All'; noise_levels = [1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]; jumpDistance = 180; species = 'dros';   for i = [1:trials] disp(sprintf('Iteration = %d\n', i)); j=2; stim_j = [j]; stim_strength_fact=1; collect_stats_long_run(species,   noise_on_synapse_type,   noise_levels, 1,        sprintf('Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_%s_jump_spiking_0_%ddeg_synapserand_%s_trials%d.mat', species, jumpDistance, noise_on_synapse_type, trials),                 'init_and_jump_bump', stim_j, stim_strength_fact, 4.0, 'spiking', jumpDistance); end
%%    trials = 10; noise_on_synapse_type = 'All'; noise_levels = [1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]; jumpDistance = 180; species = 'locust'; for i = [1:trials] disp(sprintf('Iteration = %d\n', i)); j=2; stim_j = [j]; stim_strength_fact=1; collect_stats_long_run(species,   noise_on_synapse_type,   noise_levels, 1,        sprintf('Drift_Analysis/Data/Stats/Long_Run_Sim/collect_stats_long_run_data_%s_jump_spiking_0_%ddeg_synapserand_%s_trials%d.mat', species, jumpDistance, noise_on_synapse_type, trials),                 'init_and_jump_bump', stim_j, stim_strength_fact, 4.0, 'spiking', jumpDistance); end

function collect_stats_long_run(exp_model, changed_synapse_type, percent_noise_values, display_plots, data_filename, stim_type, stim_neuron_list, stim_strength_fact, simulationDuration, spiking_or_rate_based, jumpDistance, SimulationParamsStruct)

    initial_heading_index = 2; % Angular position of the first of two pulses
    if ~exist('jumpDistance', 'var')
        jumpDistance = 180;        % Angular amount of stimulus change between first and second pulse
    end
    
    if ~exist('SimulationParamsStruct', 'var')
        SimulationParamsStruct = [];
    end
    
    % A cell array to store collected data
    data_store = {};
    
    if length(percent_noise_values) == 2 && all(percent_noise_values == [-1, -1])
        values = [10, 20, 30, 40, 50, 60, 70, 80, 100];
        values = [-values 0 values];
    else
        values = percent_noise_values;
    end
    fprintf('Noise percentage values to test:');
    fprintf(' %.1f%%', values);
    fprintf('\n\n');
    
    if ~exist('stim_type', 'var')
        stim_type = 'init_only';
    end
    
    if ~exist('stim_neuron_list', 'var')
        stim_neuron_list = [];
    end
    
    if ~exist('simulationDuration', 'var')
        simulationDuration = 60.0; % in sec
    end
    
    if ~exist('spiking_or_rate_based', 'var')
        spiking_or_rate_based = 'spiking';
    end
    
    % Stimulus specification
    % First stimulus specification
    keySet1 = {'pulseStartTime', 'pulseDuration', 'stim_neuron_list', 'maxSpikeRate', 'maxSpikeRateFactor'};
    valueSet1 = {0, 0.750, stim_neuron_list, 1, 1};
    selectedColumnsBool = containers.Map(keySet1, valueSet1);

    % Jumping bump, we only specify the time of the second pulse
    % This is a dummy map that is not actually being used. 
    keySet_jumpBarBool = {'pulseStartTime', 'pulseDuration', 'stim_neuron_list', 'jumpDistance', 'maxSpikeRate', 'maxSpikeRateFactor'};
    valueSet_jumpBarBool = {0, 0, stim_neuron_list, jumpDistance, 1, 1};
    jumpBarBool = containers.Map(keySet_jumpBarBool, valueSet_jumpBarBool);
    
    % In the case of vonMisesFixedBool the stim_neuron_list must contain 
    % one real value. Eg can be 1, 1.2, 1.5, 2, etc.
    keySet1vm = {'pulseStartTime', 'pulseDuration', 'stim_neuron_list', 'maxSpikeRate', 'maxSpikeRateFactor', 'kappa'};
    valueSet1vm = {0, 0.750, stim_neuron_list, 1, 1, pi*7/6}; % kappa=pi*7/6 results to Gaussian with 91deg FWHM
    valueSet1vm = {0, 0.750, stim_neuron_list, 1, 1, pi*2}; % kappa=pi*2 results to Gaussian with 61deg FWHM
    vonMisesFixedBool   = containers.Map(keySet1vm, valueSet1vm);

    
    % Real fruit fly flight trajectory stimulus
    if strcmp(stim_type, 'real_trajectoy_vonMises')
        % Load the fruit fly trajectory time series
        load('Fruit-fly-trajectories/fruit_fly_trajectories.mat');
        stim_neuron_center_idx_ts = trajectory_integral_headings; % Heading time series
        % Shift the headings so that 0 is moved to 180deg
        stim_neuron_center_idx_ts = mod(stim_neuron_center_idx_ts + 3, 8)+1;
        timeSeriesStartTime = min(trajectory_t); % It is 0[s] (trajectory start time)
        timeSeriesDuration  = max(trajectory_t) - min(trajectory_t); % It is 17[s] (trajectory end time)
        % Set the submap values. The same variables are used but it is not a
        % pulse then, it is a continues bump with moving heading. 
        % pulseStartTime   : trajectory start time
        % pulseDuration    : trajectory end time
        % stim_neuron_list : the heading time series (value 1-8 per time step)
        keySet1rtvm = {'pulseStartTime', 'pulseDuration', 'stim_neuron_list', 'maxSpikeRate', 'maxSpikeRateFactor', 'kappa'};
        valueSet1rtvm = {timeSeriesStartTime, timeSeriesDuration, stim_neuron_center_idx_ts, 1, 1, pi*2}; % kappa=pi*2 results to Gaussian with 61deg FWHM
        vonMisesTrajectoryBool   = containers.Map(keySet1rtvm, valueSet1rtvm);
    end
    
    if strcmp(stim_type, 'init_only')
        % Apply only an bump initialisation pulse
        keySet = {'simulationDuration', 'stimulatedNeuronGroup', 'selectedColumnsBool'};
        valueSet = {simulationDuration, 'E-PG', selectedColumnsBool};
        StimulusParamsMap = containers.Map(keySet, valueSet);
    elseif strcmp(stim_type, 'init_and_jump_bump')
        % An initialisation bump pulse and later a second pulse a specific
        % angle away from the first. Only the second pulse time is
        % specified. 
        % We use the old way of specifying the stimulus so no need to
        % specify a map here. This map is not actually used to pass values 
        % it is only used to get through the function calls that expect it.
        keySet = {'simulationDuration', 'stimulatedNeuronGroup', 'jumpBarBool'};
        valueSet = {simulationDuration, 'E-PG', jumpBarBool};
        StimulusParamsMap = containers.Map(keySet, valueSet);
    elseif strcmp(stim_type, 'init_and_unihemispheric')
        % Second stimulus specification
        keySet2 = {'pulseStartTime', 'pulseDuration', 'stim_neuron_list', 'maxSpikeRate', 'maxSpikeRateFactor'};
        valueSet2 = {5, 55, [1:8], 1, stim_strength_fact};
        ImbalanceRightBool = containers.Map(keySet2, valueSet2);
        % Apply a bump initialisation pulse and then uni-hemispheric stimulation
        keySet = {'simulationDuration', 'stimulatedNeuronGroup', 'selectedColumnsBool', 'ImbalanceRightBool'};
        valueSet = {simulationDuration, 'P-EN', selectedColumnsBool, ImbalanceRightBool};
        StimulusParamsMap = containers.Map(keySet, valueSet);
    elseif strcmp(stim_type, 'init_only_vonMises')
        % Apply only an bump initialisation pulse but having the bump width
        % properties
        keySet = {'simulationDuration', 'stimulatedNeuronGroup', 'vonMisesFixedBool'};
        valueSet = {simulationDuration, 'E-PG', vonMisesFixedBool};
        StimulusParamsMap = containers.Map(keySet, valueSet);            
    elseif strcmp(stim_type, 'real_trajectoy_vonMises')
        % Apply continusly a bump that moves as a real fruit fly heading 
        % trajectory and having the bump width properties
        keySet = {'simulationDuration', 'stimulatedNeuronGroup', 'vonMisesTrajectoryBool'};
        valueSet = {simulationDuration, 'E-PG', vonMisesTrajectoryBool};
        StimulusParamsMap = containers.Map(keySet, valueSet);            
    end
    
    for i=values
        disp(sprintf('Running simulation with noise %d%%', i));
        % Run a simulation
        if strcmp(exp_model, 'dros')
            weight_change_percent = i;
            out=dros_with_noise(changed_synapse_type, weight_change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct);
        elseif strcmp(exp_model, 'locust')
            weight_change_percent = i;
            out=locust_with_noise(changed_synapse_type, weight_change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct);
        elseif strcmp(exp_model, 'mixed_species')
            % Connectivity of drosophila with inhibitory connections
            % replaced with those of the locust model
            weight_change_percent = i;
            out=mixed_species_dros_with_noise(changed_synapse_type, weight_change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct);
        elseif strcmp(exp_model, 'dros_without_P_EG')
            weight_change_percent = i;
            out=dros_without_P_EG_with_noise(changed_synapse_type, weight_change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct);
        elseif strcmp(exp_model, 'locust_without_P_EG')
            weight_change_percent = i;
            out=locust_without_P_EG_with_noise(changed_synapse_type, weight_change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct);
        end
        
        if display_plots
            figure(); plot(out.theta_ts);
            %figure(); plot(out.r_ts);
        end
        
        % Construct a data point record
        data_record{1} = exp_model;
        data_record{2} = sparse(logical(out.out.spikeSequence));
        data_record{3} = out.n_P_EN;
        data_record{4} = out.n_P_EG;
        data_record{5} = out.n_E_PG;
        data_record{6} = out.n_Pintr;
        data_record{7} = out.con_matrix_filename;
        data_record{8} = out.con_matrix_parameter_set;
        data_record{9} = out.inhibition_distr_type;
        data_record{10} = out.inhibition_width_sigma;
        data_record{11} = out.out.P.dt;
        data_record{12} = out.simulationTime;    % Simulation duration in sec
        data_record{13} = out.changed_synapse_type;
        data_record{14} = out.change_percent;
        data_record{15} = out.con_matrix;
        % data_record{16} = out.stim_neuron_list;  % Stimulus applied to neuron IDs within the out.inputList
        data_record{16} = StimulusParamsMap('stimulatedNeuronGroup');
        data_record{17} = out.inputList;           % The input neuron group is common for all stimuli
        
        % Append stimulus entries to the cell array
        s_i = 17; % Cell array index base
        for key_cell = keys(StimulusParamsMap)
            key = key_cell{1};
            if ~strcmp(key, 'stimulatedNeuronGroup') & ~strcmp(key, 'simulationDuration')
                subMap = StimulusParamsMap(key);
                s_i=s_i+1; data_record{s_i} = key;                                     % Type of stimulus, one of 'selectedColumnsBool', 'ImbalanceRightBool'
                s_i=s_i+1; data_record{s_i} = length(subMap.keys());        % Num of entries for this key
                s_i=s_i+1; data_record{s_i} = subMap('pulseStartTime');     % Start this stimulus at this time in sec
                s_i=s_i+1; data_record{s_i} = subMap('pulseDuration');      % This stimulus duration in sec
                s_i=s_i+1; data_record{s_i} = subMap('stim_neuron_list');   % Stimulus applied to these neuron IDs within the out.inputList
                s_i=s_i+1; data_record{s_i} = subMap('maxSpikeRate');       % maximum/peak spike rate for this stimulus
                s_i=s_i+1; data_record{s_i} = subMap('maxSpikeRateFactor'); % Multiplying factor specifying portion of the maxSpikeRate that should be applied
                if isKey(subMap, 'jumpDistance')
                    s_i=s_i+1; data_record{s_i} = subMap('jumpDistance');   % Multiplying factor specifying portion of the maxSpikeRate that should be applied
                end
            end
        end
        
        s_i=s_i+1; data_record{s_i} = out.theta_ts; % The azimuth of the bump centre time series
        s_i=s_i+1; data_record{s_i} = out.r_ts;     % The length of the vector representing the bump time series
        
        % Store the record
        data_store{1} = data_record;
        add_data_to_file(data_filename, data_store);
    end
end




%% Adds new data entry to the file
function add_data_to_file(data_filename, data)
    if exist(data_filename, 'file') == 2
        % Loads a variable cell_array which is a cell array
        load(data_filename);
    else
        cell_array = [];
    end
    
    for i=1:size(data, 2)
        cell_array = [cell_array; data{i}];
    end

    save(data_filename, 'cell_array', '-v7.3');
end

%% Drosophila connectivity With P-EG neurons
function out=dros_with_noise(changed_synapse_type, change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct)
    fixedThreshold = false;
    jumpDistance = jumpDistance;
    currentStrenth = 0.0;
    pulseLengthMulti = 3;

    maxSpikeRate = 170; % 170 Hz is the original value used.

    % Drosophila corrected (2018) connectivity matrix including P-EG neurons
    con_matrix_filename = 'connectivity_matrix_drosophila_mine_case_5_9cols_labels1.mat';
    load(con_matrix_filename);
    disp('Drosophila corrected (2018) connectivity matrix including P-EG neurons')
    
    if strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'E-PG')
        %           [       35      :       52           ];
        inputList = [n_P_EN+n_P_EG+1:n_P_EN+n_P_EG+n_E_PG]; % E-PG neurons indexes for Drosophila with 9 columns connected
    elseif strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'P-EN')
        %           [        1      :       16           ];
        inputList = [1:n_P_EN]; % P-EN neurons indexes for Drosophila with 9 columns connected
    end
    
    % Find matrix elements with these label values
    plus_1s_ind    = con_matrix == 1;plus_3s_ind    = con_matrix == 3;plus_1_1s_ind  = con_matrix == 1.1;
    plus_0_9s_ind = con_matrix == 0.9;
    plus_0_5s_ind  = con_matrix == 0.5;
    plus_0_1s_ind = con_matrix == 0.1;
    minus_1s_ind   = con_matrix == -1;minus_1_1s_ind = con_matrix == -1.1;minus_0_5s_ind  = con_matrix == -0.5;
    minus_0_9s_ind = con_matrix == -0.9;
    minus_0_1s_ind = con_matrix == -0.1;

    % Put ones on the P_EN output synaptic value indexes
    base_n_x = 0;       % Start index
    range_n_x = n_P_EN; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EN_E_PG_output_synapses_indx,P_EN_E_PG_output_synapses_indx_left,P_EN_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);
    
    base_n_x = n_P_EN;  % Start index
    range_n_x = n_P_EG; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EG_E_PG_output_synapses_indx,P_EG_E_PG_output_synapses_indx_left,P_EG_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = 0;               % Start index
    range_n_y = n_P_EN;         % How many indeces from start
    [E_PG_P_EN_output_synapses_indx,E_PG_P_EN_output_synapses_indx_left,E_PG_P_EN_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = n_P_EN;          % Start index
    range_n_y = n_P_EG;         % How many indeces from start
    [E_PG_P_EG_output_synapses_indx,E_PG_P_EG_output_synapses_indx_left,E_PG_P_EG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG;  % Start index
    range_n_x = n_E_PG;          % How many indeces from start
    base_n_y = n_P_EN + n_P_EG + n_E_PG; % Start index
    range_n_y = n_Pintr;                 % How many indeces from start
    [E_PG_Pintr_output_synapses_indx,E_PG_Pintr_output_synapses_indx_left,E_PG_Pintr_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);
    
    base_n_x = n_P_EN + n_P_EG + n_E_PG; % Start index
    range_n_x = n_Pintr;                 % How many indeces from start
    base_n_y = 0;                   % Start index
    range_n_y = n_P_EN;                  % How many indeces from start
    [Pintr_P_EN_output_synapses_indx,Pintr_P_EN_output_synapses_indx_left,Pintr_P_EN_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG + n_E_PG; % Start index
    range_n_x = n_Pintr;                 % How many indeces from start
    base_n_y = n_P_EN;                   % Start index
    range_n_y = n_P_EG;                  % How many indeces from start
    [Pintr_P_EG_output_synapses_indx,Pintr_P_EG_output_synapses_indx_left,Pintr_P_EG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG + n_E_PG; % Start index
    range_n_x = n_Pintr;                 % How many indeces from start
    base_n_y = n_P_EN + n_P_EG + n_E_PG;                   % Start index
    range_n_y = n_Pintr;                  % How many indeces from start
    [Pintr_Pintr_output_synapses_indx,Pintr_Pintr_output_synapses_indx_left,Pintr_Pintr_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);
    
    
    % <entry version> <fitness> <bumpLocPre> <bumpLocPost> <bumpWidthPre> <final bump width> <E-PG peak spike rate> <5 or 8 x weight parameters>
    inhibition_distr_type = 'uniform';
    inhibition_width_sigma = 0.0;
    parameter_set_drosophila = [5	0.0234854383940386	1.50000000000000	6	0.785398163397448	1.57079632679490	190	6.07203896674184	47.9173822984772	-35.6772759590409	-19.3358579258735	99.9680554605268];

    if size(parameter_set_drosophila, 2) == 5
        [plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
    elseif size(parameter_set_drosophila, 2) == 6
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        plus_3s_value    = plus_1s_value * 3;
    elseif size(parameter_set_drosophila, 2) == 7
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
    elseif size(parameter_set_drosophila, 2) == 9
        [version,fitness,bump_width,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        maxSpikeRate = peak_spike_rate_E_PG;
    elseif size(parameter_set_drosophila, 2) == 12
        [version,fitness,bump_loc_pre,bump_loc_post,bump_width_pre,bump_width_post,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        bump_width = bump_width_post;
        maxSpikeRate = peak_spike_rate_E_PG;
    end

    % If the connectivity matrix contains -0.9 values then is sinusoidal
    % otherwise is the hardwired ratio
    if any(any(minus_0_9s_ind))
        half_mult = 0.5;
    else
        half_mult = 0.707;
    end

    minus_0_9s_value = minus_1_1s_value * 0.85355339;
    minus_0_5s_value = minus_1_1s_value * half_mult;
    minus_0_1s_value = minus_1_1s_value * 0.14644661;

    plus_0_9s_value = plus_1_1s_value * 0.85355339;
    plus_0_5s_value = plus_1_1s_value * half_mult;
    plus_0_1s_value = plus_1_1s_value * 0.14644661;

    % Print the parameter values
    [plus_1s_value plus_1_1s_value plus_0_5s_value minus_1s_value minus_1_1s_value minus_0_5s_value plus_3s_value];
    con_matrix(plus_1s_ind)    = plus_1s_value;con_matrix(plus_3s_ind)    = plus_3s_value;con_matrix(plus_1_1s_ind)  = plus_1_1s_value;con_matrix(plus_0_9s_ind)  = plus_0_9s_value;con_matrix(plus_0_5s_ind)  = plus_0_5s_value;con_matrix(plus_0_1s_ind)  = plus_0_1s_value;con_matrix(minus_1s_ind)   = minus_1s_value;con_matrix(minus_1_1s_ind) = minus_1_1s_value;con_matrix(minus_0_9s_ind) = minus_0_9s_value;con_matrix(minus_0_5s_ind) = minus_0_5s_value;con_matrix(minus_0_1s_ind) = minus_0_1s_value;

    % Add noise to synaptic weights
    % Alter synaptic weights of all synapses
    if strcmp(changed_synapse_type, 'All')
        % Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix)) .* con_matrix * (change_percent / 100);
        con_matrix = con_matrix + added_noise;
    end
    % Alter synaptic weights of one type of synapses on one side or randomly on both sides
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_left')
        con_matrix(P_EN_E_PG_output_synapses_indx_left) = con_matrix(P_EN_E_PG_output_synapses_indx_left) + con_matrix(P_EN_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_right')
        con_matrix(P_EN_E_PG_output_synapses_indx_right) = con_matrix(P_EN_E_PG_output_synapses_indx_right) + con_matrix(P_EN_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EN_E_PG_output_synapses_indx))) .* con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_alleq')
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_left')
        con_matrix(P_EG_E_PG_output_synapses_indx_left) = con_matrix(P_EG_E_PG_output_synapses_indx_left) + con_matrix(P_EG_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_right')
        con_matrix(P_EG_E_PG_output_synapses_indx_right) = con_matrix(P_EG_E_PG_output_synapses_indx_right) + con_matrix(P_EG_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EG_E_PG_output_synapses_indx))) .* con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_alleq')
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_left')
        con_matrix(E_PG_P_EN_output_synapses_indx_left) = con_matrix(E_PG_P_EN_output_synapses_indx_left) + con_matrix(E_PG_P_EN_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_right')
        con_matrix(E_PG_P_EN_output_synapses_indx_right) = con_matrix(E_PG_P_EN_output_synapses_indx_right) + con_matrix(E_PG_P_EN_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EN_output_synapses_indx))) .* con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_alleq')
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_left')
        con_matrix(E_PG_P_EG_output_synapses_indx_left) = con_matrix(E_PG_P_EG_output_synapses_indx_left) + con_matrix(E_PG_P_EG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_right')
        con_matrix(E_PG_P_EG_output_synapses_indx_right) = con_matrix(E_PG_P_EG_output_synapses_indx_right) + con_matrix(E_PG_P_EG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EG_output_synapses_indx))) .* con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_alleq')
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_left')
        con_matrix(E_PG_Pintr_output_synapses_indx_left) = con_matrix(E_PG_Pintr_output_synapses_indx_left) + con_matrix(E_PG_Pintr_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_right')
        con_matrix(E_PG_Pintr_output_synapses_indx_right) = con_matrix(E_PG_Pintr_output_synapses_indx_right) + con_matrix(E_PG_Pintr_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_Pintr_output_synapses_indx))) .* con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_alleq')
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
    end
        
    if strcmp(changed_synapse_type, 'Pintr_to_P_EN_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(Pintr_P_EN_output_synapses_indx))) .* con_matrix(Pintr_P_EN_output_synapses_indx) * (change_percent / 100);
        con_matrix(Pintr_P_EN_output_synapses_indx) = con_matrix(Pintr_P_EN_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'Pintr_to_P_EN_alleq')
        con_matrix(Pintr_P_EN_output_synapses_indx) = con_matrix(Pintr_P_EN_output_synapses_indx) + con_matrix(Pintr_P_EN_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'Pintr_to_P_EG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(Pintr_P_EG_output_synapses_indx))) .* con_matrix(Pintr_P_EG_output_synapses_indx) * (change_percent / 100);
        con_matrix(Pintr_P_EG_output_synapses_indx) = con_matrix(Pintr_P_EG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'Pintr_to_P_EG_alleq')
        con_matrix(Pintr_P_EG_output_synapses_indx) = con_matrix(Pintr_P_EG_output_synapses_indx) + con_matrix(Pintr_P_EG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'Pintr_to_Pintr_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(Pintr_Pintr_output_synapses_indx))) .* con_matrix(Pintr_Pintr_output_synapses_indx) * (change_percent / 100);
        con_matrix(Pintr_Pintr_output_synapses_indx) = con_matrix(Pintr_Pintr_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'Pintr_to_Pintr_alleq')
        con_matrix(Pintr_Pintr_output_synapses_indx) = con_matrix(Pintr_Pintr_output_synapses_indx) + con_matrix(Pintr_Pintr_output_synapses_indx) * (change_percent / 100);
    end
    
    % Update the maxSpikeRate values with the ones read from the file
    if isKey(StimulusParamsMap, 'selectedColumnsBool')
        paramsMap = StimulusParamsMap('selectedColumnsBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'jumpBarBool')
        paramsMap = StimulusParamsMap('jumpBarBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'ImbalanceRightBool')
        paramsMap = StimulusParamsMap('ImbalanceRightBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesFixedBool')
        paramsMap = StimulusParamsMap('vonMisesFixedBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesTrajectoryBool')
        paramsMap = StimulusParamsMap('vonMisesTrajectoryBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    
    % Run the network simulation and get bump width at the end
    bump = runNetworkSimulation(con_matrix, display_plots, fixedThreshold, currentStrenth, jumpDistance, pulseLengthMulti, maxSpikeRate, inputList, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, SimulationParamsStruct);
    bumpWidth = bump.bumpWidth;
    % Print the bump width
    %bumpWidth * 180 / pi

    %bump.bumpLocPre
    %bump.bumpLocPost

    %bump.bumpWidthPre  * 180 / pi
    %bump.bumpWidthPost * 180 / pi

    %bump.spikeRate % Display the spike rates of the neurons

    out = bump;
    out.n_P_EN  = n_P_EN;
    out.n_P_EG  = n_P_EG;
    out.n_E_PG  = n_E_PG;
    out.n_Pintr = n_Pintr;
    out.con_matrix_filename = con_matrix_filename;
    out.con_matrix_parameter_set = parameter_set_drosophila;
    out.inhibition_distr_type = inhibition_distr_type;
    out.inhibition_width_sigma = inhibition_width_sigma;
    out.changed_synapse_type = changed_synapse_type;
    out.change_percent = change_percent;
    out.con_matrix = con_matrix;
    out.inputList = inputList;
    out.maxSpikeRate = maxSpikeRate;
    
end


%% Locust connectivity With P-EG neurons
function out=locust_with_noise(changed_synapse_type, change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct)
    fixedThreshold = false;
    jumpDistance = jumpDistance;
    currentStrenth = 0.0;
    pulseLengthMulti = 3;
    maxSpikeRate = 170; % 170 Hz is the original value used.

    % Locust connectivity matrix including P-EG neurons
    con_matrix_filename = 'connectivity_matrix_locust_mine_case_6_with3s_labels_1.mat';
    load(con_matrix_filename);
    disp('Locust connectivity matrix including P-EG neurons')

    if strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'E-PG')
        %           [       33      :       48           ];
        inputList = [n_P_EN+n_P_EG+1:n_P_EN+n_P_EG+n_E_PG]; % E-PG neurons indexes for locust with 8 columns connected
    elseif strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'P-EN')
        %           [        1      :       16           ];
        inputList = [1:n_P_EN]; % P-EN neurons indexes for locust with 8 columns connected
    end

    % Find matrix elements with these label values
    plus_1s_ind    = con_matrix == 1;plus_3s_ind    = con_matrix == 3;plus_1_1s_ind  = con_matrix == 1.1;
    plus_0_9s_ind = con_matrix == 0.9;
    plus_0_5s_ind  = con_matrix == 0.5;
    plus_0_1s_ind = con_matrix == 0.1;
    minus_1s_ind   = con_matrix == -1;minus_1_1s_ind = con_matrix == -1.1;minus_0_5s_ind  = con_matrix == -0.5;
    minus_0_9s_ind = con_matrix == -0.9;
    minus_0_1s_ind = con_matrix == -0.1;

    % Put ones on the P_EN output synaptic value indexes
    base_n_x = 0;       % Start index
    range_n_x = n_P_EN; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EN_E_PG_output_synapses_indx,P_EN_E_PG_output_synapses_indx_left,P_EN_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);
    
    base_n_x = n_P_EN;  % Start index
    range_n_x = n_P_EG; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EG_E_PG_output_synapses_indx,P_EG_E_PG_output_synapses_indx_left,P_EG_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = 0;               % Start index
    range_n_y = n_P_EN;         % How many indeces from start
    [E_PG_P_EN_output_synapses_indx,E_PG_P_EN_output_synapses_indx_left,E_PG_P_EN_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = n_P_EN;          % Start index
    range_n_y = n_P_EG;         % How many indeces from start
    [E_PG_P_EG_output_synapses_indx,E_PG_P_EG_output_synapses_indx_left,E_PG_P_EG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG;  % Start index
    range_n_x = n_E_PG;          % How many indeces from start
    base_n_y = n_P_EN + n_P_EG + n_E_PG; % Start index
    range_n_y = n_Pintr;                 % How many indeces from start
    [E_PG_Pintr_output_synapses_indx,E_PG_Pintr_output_synapses_indx_left,E_PG_Pintr_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    % <entry version> <fitness> <bumpLocPre> <bumpLocPost> <bumpWidthPre> <final bump width> <E-PG peak spike rate> <5 or 8 x weight parameters>
    inhibition_distr_type = 'gaussian';
    inhibition_width_sigma = 0.4;
    parameter_set_locust = [  5	0.0355556337417948	2	5	1.57079632679490	1.57079632679490	180	10.8341743712825	30.7047329203947	-42.1014161760853	-9.78202338418053	20.2648078601249];
    
    if size(parameter_set_locust, 2) == 5
        [plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
    elseif size(parameter_set_locust, 2) == 6
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
        plus_3s_value    = plus_1s_value * 3;
    elseif size(parameter_set_locust, 2) == 7
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
    elseif size(parameter_set_locust, 2) == 9
        [version,fitness,bump_width,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
        maxSpikeRate = peak_spike_rate_E_PG;
    elseif size(parameter_set_locust, 2) == 12
        [version,fitness,bump_loc_pre,bump_loc_post,bump_width_pre,bump_width_post,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
        bump_width = bump_width_post;
        maxSpikeRate = peak_spike_rate_E_PG;
    end

    num_of_samples = n_E_PG / 2;
    if num_of_samples == 8
        num_of_samples = 9; % That is because with 8 samples it does not approximate Tom's values
    end
    inhibition_dist = get_inhib_distr_vector(inhibition_distr_type, num_of_samples, inhibition_width_sigma);
    
    % If the connectivity matrix contains -0.9 values then is sinusoidal
    % otherwise is the hardwired ratio
    if any(any(minus_0_9s_ind))
        half_mult = inhibition_dist(3); % 0.5;
    else
        half_mult = 0.707;
    end

    minus_0_9s_value = minus_1_1s_value *  inhibition_dist(4);
    minus_0_5s_value = minus_1_1s_value * half_mult;
    minus_0_1s_value = minus_1_1s_value * inhibition_dist(2);

    plus_0_9s_value = plus_1_1s_value * inhibition_dist(4);
    plus_0_5s_value = plus_1_1s_value * half_mult;
    plus_0_1s_value = plus_1_1s_value * inhibition_dist(2);

    % Print the parameter values
    [plus_1s_value plus_1_1s_value plus_0_5s_value minus_1s_value minus_1_1s_value minus_0_5s_value plus_3s_value];
    con_matrix(plus_1s_ind)    = plus_1s_value;con_matrix(plus_3s_ind)    = plus_3s_value;con_matrix(plus_1_1s_ind)  = plus_1_1s_value;con_matrix(plus_0_9s_ind)  = plus_0_9s_value;con_matrix(plus_0_5s_ind)  = plus_0_5s_value;con_matrix(plus_0_1s_ind)  = plus_0_1s_value;con_matrix(minus_1s_ind)   = minus_1s_value;con_matrix(minus_1_1s_ind) = minus_1_1s_value;con_matrix(minus_0_9s_ind) = minus_0_9s_value;con_matrix(minus_0_5s_ind) = minus_0_5s_value;con_matrix(minus_0_1s_ind) = minus_0_1s_value;

    % Add noise to synaptic weights
    % Alter synaptic weights of all synapses
    if strcmp(changed_synapse_type, 'All')
        % Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix)) .* con_matrix * (change_percent / 100);
        con_matrix = con_matrix + added_noise;
    end
    % Alter synaptic weights of one type of synapses on one side or randomly on both sides
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_left')
        con_matrix(P_EN_E_PG_output_synapses_indx_left) = con_matrix(P_EN_E_PG_output_synapses_indx_left) + con_matrix(P_EN_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_right')
        con_matrix(P_EN_E_PG_output_synapses_indx_right) = con_matrix(P_EN_E_PG_output_synapses_indx_right) + con_matrix(P_EN_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EN_E_PG_output_synapses_indx))) .* con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_alleq')
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_left')
        con_matrix(P_EG_E_PG_output_synapses_indx_left) = con_matrix(P_EG_E_PG_output_synapses_indx_left) + con_matrix(P_EG_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_right')
        con_matrix(P_EG_E_PG_output_synapses_indx_right) = con_matrix(P_EG_E_PG_output_synapses_indx_right) + con_matrix(P_EG_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EG_E_PG_output_synapses_indx))) .* con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_alleq')
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_left')
        con_matrix(E_PG_P_EN_output_synapses_indx_left) = con_matrix(E_PG_P_EN_output_synapses_indx_left) + con_matrix(E_PG_P_EN_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_right')
        con_matrix(E_PG_P_EN_output_synapses_indx_right) = con_matrix(E_PG_P_EN_output_synapses_indx_right) + con_matrix(E_PG_P_EN_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EN_output_synapses_indx))) .* con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_alleq')
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_left')
        con_matrix(E_PG_P_EG_output_synapses_indx_left) = con_matrix(E_PG_P_EG_output_synapses_indx_left) + con_matrix(E_PG_P_EG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_right')
        con_matrix(E_PG_P_EG_output_synapses_indx_right) = con_matrix(E_PG_P_EG_output_synapses_indx_right) + con_matrix(E_PG_P_EG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EG_output_synapses_indx))) .* con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_alleq')
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_left')
        con_matrix(E_PG_Pintr_output_synapses_indx_left) = con_matrix(E_PG_Pintr_output_synapses_indx_left) + con_matrix(E_PG_Pintr_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_right')
        con_matrix(E_PG_Pintr_output_synapses_indx_right) = con_matrix(E_PG_Pintr_output_synapses_indx_right) + con_matrix(E_PG_Pintr_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_Pintr_output_synapses_indx))) .* con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_alleq')
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
    end
    
    % Update the maxSpikeRate values with the ones read from the file
    if isKey(StimulusParamsMap, 'selectedColumnsBool')
        paramsMap = StimulusParamsMap('selectedColumnsBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'jumpBarBool')
        paramsMap = StimulusParamsMap('jumpBarBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'ImbalanceRightBool')
        paramsMap = StimulusParamsMap('ImbalanceRightBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesFixedBool')
        paramsMap = StimulusParamsMap('vonMisesFixedBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesTrajectoryBool')
        paramsMap = StimulusParamsMap('vonMisesTrajectoryBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    
    % Run the network simulation and get bump width at the end
    bump = runNetworkSimulation(con_matrix, display_plots, fixedThreshold, currentStrenth, jumpDistance, pulseLengthMulti, maxSpikeRate, inputList, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, SimulationParamsStruct);
    bumpWidth = bump.bumpWidth;
    % Print the bump width
    % bumpWidth * 180 / pi

    % bump.bumpLocPre
    % bump.bumpLocPost

    % bump.bumpWidthPre  * 180 / pi
    % bump.bumpWidthPost * 180 / pi

    % bump.spikeRate % Display the spike rates of the neurons
    
    out = bump;
    out.n_P_EN  = n_P_EN;
    out.n_P_EG  = n_P_EG;
    out.n_E_PG  = n_E_PG;
    out.n_Pintr = n_Pintr;
    out.con_matrix_filename = con_matrix_filename;
    out.con_matrix_parameter_set = parameter_set_locust;
    out.inhibition_distr_type = inhibition_distr_type;
    out.inhibition_width_sigma = inhibition_width_sigma;
    out.changed_synapse_type = changed_synapse_type;
    out.change_percent = change_percent;
    out.con_matrix = con_matrix;
    out.inputList = inputList;
    out.maxSpikeRate = maxSpikeRate;
    
end


%% Drosophila connectivity Without P-EG neurons
function out=dros_without_P_EG_with_noise(changed_synapse_type, change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct)
    fixedThreshold = false;
    jumpDistance = jumpDistance;
    currentStrenth = 0.0;
    pulseLengthMulti = 3;
    maxSpikeRate = 170; % 170 Hz is the original value used.

    % Drosophila corrected (2018) connectivity matrix without P-EG neurons
    con_matrix_filename = 'connectivity_matrix_drosophila_mine_case_5_9cols_withoutP_EG_labels1.mat';
    load(con_matrix_filename);
    disp('Drosophila corrected (2018) connectivity matrix without P-EG neurons')
    
    if strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'E-PG')
        %           [       35      :       52           ];
        inputList = [n_P_EN+n_P_EG+1:n_P_EN+n_P_EG+n_E_PG]; % E-PG neurons indexes for Drosophila with 9 columns connected
    elseif strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'P-EN')
        %           [        1      :       16           ];
        inputList = [1:n_P_EN]; % P-EN neurons indexes for Drosophila with 9 columns connected
    end
    
    % Find matrix elements with these label values
    plus_1s_ind    = con_matrix == 1;plus_3s_ind    = con_matrix == 3;plus_1_1s_ind  = con_matrix == 1.1;
    plus_0_9s_ind = con_matrix == 0.9;
    plus_0_5s_ind  = con_matrix == 0.5;
    plus_0_1s_ind = con_matrix == 0.1;
    minus_1s_ind   = con_matrix == -1;minus_1_1s_ind = con_matrix == -1.1;minus_0_5s_ind  = con_matrix == -0.5;
    minus_0_9s_ind = con_matrix == -0.9;
    minus_0_1s_ind = con_matrix == -0.1;

    % Put ones on the P_EN output synaptic value indexes
    base_n_x = 0;       % Start index
    range_n_x = n_P_EN; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EN_E_PG_output_synapses_indx,P_EN_E_PG_output_synapses_indx_left,P_EN_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);
    
    base_n_x = n_P_EN;  % Start index
    range_n_x = n_P_EG; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EG_E_PG_output_synapses_indx,P_EG_E_PG_output_synapses_indx_left,P_EG_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = 0;               % Start index
    range_n_y = n_P_EN;         % How many indeces from start
    [E_PG_P_EN_output_synapses_indx,E_PG_P_EN_output_synapses_indx_left,E_PG_P_EN_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = n_P_EN;          % Start index
    range_n_y = n_P_EG;         % How many indeces from start
    [E_PG_P_EG_output_synapses_indx,E_PG_P_EG_output_synapses_indx_left,E_PG_P_EG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG;  % Start index
    range_n_x = n_E_PG;          % How many indeces from start
    base_n_y = n_P_EN + n_P_EG + n_E_PG; % Start index
    range_n_y = n_Pintr;                 % How many indeces from start
    [E_PG_Pintr_output_synapses_indx,E_PG_Pintr_output_synapses_indx_left,E_PG_Pintr_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    % <entry version> <fitness> <bumpLocPre> <bumpLocPost> <bumpWidthPre> <final bump width> <E-PG peak spike rate> <5 or 8 x weight parameters>
    inhibition_distr_type = 'uniform';
    inhibition_width_sigma = 0.0;
    parameter_set_drosophila = [5	0.0242425305349241	2	7	1.57079632679490	1.57079632679490	210	58.6429952973772	77.2366215106736	-41.8245611289463	-67.5674037322441	64.5831383390314]; % 

    if size(parameter_set_drosophila, 2) == 5
        [plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
    elseif size(parameter_set_drosophila, 2) == 6
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        plus_3s_value    = plus_1s_value * 3;
    elseif size(parameter_set_drosophila, 2) == 7
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
    elseif size(parameter_set_drosophila, 2) == 9
        [version,fitness,bump_width,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        maxSpikeRate = peak_spike_rate_E_PG;
    elseif size(parameter_set_drosophila, 2) == 12
        [version,fitness,bump_loc_pre,bump_loc_post,bump_width_pre,bump_width_post,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        bump_width = bump_width_post;
        maxSpikeRate = peak_spike_rate_E_PG;
    end

    % If the connectivity matrix contains -0.9 values then is sinusoidal
    % otherwise is the hardwired ratio
    if any(any(minus_0_9s_ind))
        half_mult = 0.5;
    else
        half_mult = 0.707;
    end

    minus_0_9s_value = minus_1_1s_value * 0.85355339;
    minus_0_5s_value = minus_1_1s_value * half_mult;
    minus_0_1s_value = minus_1_1s_value * 0.14644661;

    plus_0_9s_value = plus_1_1s_value * 0.85355339;
    plus_0_5s_value = plus_1_1s_value * half_mult;
    plus_0_1s_value = plus_1_1s_value * 0.14644661;

    % Print the parameter values
    [plus_1s_value plus_1_1s_value plus_0_5s_value minus_1s_value minus_1_1s_value minus_0_5s_value plus_3s_value];
    con_matrix(plus_1s_ind)    = plus_1s_value;con_matrix(plus_3s_ind)    = plus_3s_value;con_matrix(plus_1_1s_ind)  = plus_1_1s_value;con_matrix(plus_0_9s_ind)  = plus_0_9s_value;con_matrix(plus_0_5s_ind)  = plus_0_5s_value;con_matrix(plus_0_1s_ind)  = plus_0_1s_value;con_matrix(minus_1s_ind)   = minus_1s_value;con_matrix(minus_1_1s_ind) = minus_1_1s_value;con_matrix(minus_0_9s_ind) = minus_0_9s_value;con_matrix(minus_0_5s_ind) = minus_0_5s_value;con_matrix(minus_0_1s_ind) = minus_0_1s_value;

    % Add noise to synaptic weights
    % Alter synaptic weights of all synapses
    if strcmp(changed_synapse_type, 'All')
        % Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix)) .* con_matrix * (change_percent / 100);
        con_matrix = con_matrix + added_noise;
    end
    % Alter synaptic weights of one type of synapses on one side or randomly on both sides
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_left')
        con_matrix(P_EN_E_PG_output_synapses_indx_left) = con_matrix(P_EN_E_PG_output_synapses_indx_left) + con_matrix(P_EN_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_right')
        con_matrix(P_EN_E_PG_output_synapses_indx_right) = con_matrix(P_EN_E_PG_output_synapses_indx_right) + con_matrix(P_EN_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EN_E_PG_output_synapses_indx))) .* con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_alleq')
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_left')
        con_matrix(P_EG_E_PG_output_synapses_indx_left) = con_matrix(P_EG_E_PG_output_synapses_indx_left) + con_matrix(P_EG_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_right')
        con_matrix(P_EG_E_PG_output_synapses_indx_right) = con_matrix(P_EG_E_PG_output_synapses_indx_right) + con_matrix(P_EG_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EG_E_PG_output_synapses_indx))) .* con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_alleq')
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_left')
        con_matrix(E_PG_P_EN_output_synapses_indx_left) = con_matrix(E_PG_P_EN_output_synapses_indx_left) + con_matrix(E_PG_P_EN_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_right')
        con_matrix(E_PG_P_EN_output_synapses_indx_right) = con_matrix(E_PG_P_EN_output_synapses_indx_right) + con_matrix(E_PG_P_EN_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EN_output_synapses_indx))) .* con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_alleq')
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_left')
        con_matrix(E_PG_P_EG_output_synapses_indx_left) = con_matrix(E_PG_P_EG_output_synapses_indx_left) + con_matrix(E_PG_P_EG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_right')
        con_matrix(E_PG_P_EG_output_synapses_indx_right) = con_matrix(E_PG_P_EG_output_synapses_indx_right) + con_matrix(E_PG_P_EG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EG_output_synapses_indx))) .* con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_alleq')
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_left')
        con_matrix(E_PG_Pintr_output_synapses_indx_left) = con_matrix(E_PG_Pintr_output_synapses_indx_left) + con_matrix(E_PG_Pintr_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_right')
        con_matrix(E_PG_Pintr_output_synapses_indx_right) = con_matrix(E_PG_Pintr_output_synapses_indx_right) + con_matrix(E_PG_Pintr_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_Pintr_output_synapses_indx))) .* con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_alleq')
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
    end
    
    % Update the maxSpikeRate values with the ones read from the file
    if isKey(StimulusParamsMap, 'selectedColumnsBool')
        paramsMap = StimulusParamsMap('selectedColumnsBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'jumpBarBool')
        paramsMap = StimulusParamsMap('jumpBarBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'ImbalanceRightBool')
        paramsMap = StimulusParamsMap('ImbalanceRightBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesFixedBool')
        paramsMap = StimulusParamsMap('vonMisesFixedBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesTrajectoryBool')
        paramsMap = StimulusParamsMap('vonMisesTrajectoryBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    
    % Run the network simulation and get bump width at the end
    bump = runNetworkSimulation(con_matrix, display_plots, fixedThreshold, currentStrenth, jumpDistance, pulseLengthMulti, maxSpikeRate, inputList, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, SimulationParamsStruct);
    bumpWidth = bump.bumpWidth;
    % Print the bump width
    %bumpWidth * 180 / pi

    %bump.bumpLocPre
    %bump.bumpLocPost

    %bump.bumpWidthPre  * 180 / pi
    %bump.bumpWidthPost * 180 / pi

    %bump.spikeRate % Display the spike rates of the neurons

    out = bump;
    out.n_P_EN  = n_P_EN;
    out.n_P_EG  = n_P_EG;
    out.n_E_PG  = n_E_PG;
    out.n_Pintr = n_Pintr;
    out.con_matrix_filename = con_matrix_filename;
    out.con_matrix_parameter_set = parameter_set_drosophila;
    out.inhibition_distr_type = inhibition_distr_type;
    out.inhibition_width_sigma = inhibition_width_sigma;
    out.changed_synapse_type = changed_synapse_type;
    out.change_percent = change_percent;
    out.con_matrix = con_matrix;
    out.inputList = inputList;
    out.maxSpikeRate = maxSpikeRate;
    
end


%% Locust connectivity Without P-EG neurons
function out=locust_without_P_EG_with_noise(changed_synapse_type, change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct)
    fixedThreshold = false;
    jumpDistance = jumpDistance;
    currentStrenth = 0.0;
    pulseLengthMulti = 3;
    maxSpikeRate = 170; % 170 Hz is the original value used.

    % Locust connectivity matrix without P-EG neurons
    con_matrix_filename = 'connectivity_matrix_locust_mine_case_6_with3s_without_P-EG_labels_1.mat';
    load(con_matrix_filename);
    disp('Locust connectivity matrix without P-EG neurons')

    if strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'E-PG')
        %           [       33      :       48           ];
        inputList = [n_P_EN+n_P_EG+1:n_P_EN+n_P_EG+n_E_PG]; % E-PG neurons indexes for locust with 8 columns connected
    elseif strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'P-EN')
        %           [        1      :       16           ];
        inputList = [1:n_P_EN]; % P-EN neurons indexes for locust with 8 columns connected
    end

    % Find matrix elements with these label values
    plus_1s_ind    = con_matrix == 1;plus_3s_ind    = con_matrix == 3;plus_1_1s_ind  = con_matrix == 1.1;
    plus_0_9s_ind = con_matrix == 0.9;
    plus_0_5s_ind  = con_matrix == 0.5;
    plus_0_1s_ind = con_matrix == 0.1;
    minus_1s_ind   = con_matrix == -1;minus_1_1s_ind = con_matrix == -1.1;minus_0_5s_ind  = con_matrix == -0.5;
    minus_0_9s_ind = con_matrix == -0.9;
    minus_0_1s_ind = con_matrix == -0.1;

    % Put ones on the P_EN output synaptic value indexes
    base_n_x = 0;       % Start index
    range_n_x = n_P_EN; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EN_E_PG_output_synapses_indx,P_EN_E_PG_output_synapses_indx_left,P_EN_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);
    
    base_n_x = n_P_EN;  % Start index
    range_n_x = n_P_EG; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EG_E_PG_output_synapses_indx,P_EG_E_PG_output_synapses_indx_left,P_EG_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = 0;               % Start index
    range_n_y = n_P_EN;         % How many indeces from start
    [E_PG_P_EN_output_synapses_indx,E_PG_P_EN_output_synapses_indx_left,E_PG_P_EN_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = n_P_EN;          % Start index
    range_n_y = n_P_EG;         % How many indeces from start
    [E_PG_P_EG_output_synapses_indx,E_PG_P_EG_output_synapses_indx_left,E_PG_P_EG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG;  % Start index
    range_n_x = n_E_PG;          % How many indeces from start
    base_n_y = n_P_EN + n_P_EG + n_E_PG; % Start index
    range_n_y = n_Pintr;                 % How many indeces from start
    [E_PG_Pintr_output_synapses_indx,E_PG_Pintr_output_synapses_indx_left,E_PG_Pintr_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    % <entry version> <fitness> <bumpLocPre> <bumpLocPost> <bumpWidthPre> <final bump width> <E-PG peak spike rate> <5 or 8 x weight parameters>
    inhibition_distr_type = 'gaussian';
    inhibition_width_sigma = 0.4;
    parameter_set_locust = [  5	7.79478332883703e-08	2	6	1.57079632679490	1.57079632679490	200	20.0616192090804	90.3585470945102	-34.1854743836289	-17.3090658712475	49.4552854091619  ];

    if size(parameter_set_locust, 2) == 5
        [plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
    elseif size(parameter_set_locust, 2) == 6
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
        plus_3s_value    = plus_1s_value * 3;
    elseif size(parameter_set_locust, 2) == 7
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
    elseif size(parameter_set_locust, 2) == 9
        [version,fitness,bump_width,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
        maxSpikeRate = peak_spike_rate_E_PG;
    elseif size(parameter_set_locust, 2) == 12
        [version,fitness,bump_loc_pre,bump_loc_post,bump_width_pre,bump_width_post,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_locust));
        bump_width = bump_width_post;
        maxSpikeRate = peak_spike_rate_E_PG;
    end

    num_of_samples = n_E_PG / 2;
    if num_of_samples == 8
        num_of_samples = 9; % That is because with 8 samples it does not approximate Tom's values
    end
    inhibition_dist = get_inhib_distr_vector(inhibition_distr_type, num_of_samples, inhibition_width_sigma);
    
    % If the connectivity matrix contains -0.9 values then is sinusoidal
    % otherwise is the hardwired ratio
    if any(any(minus_0_9s_ind))
        half_mult = inhibition_dist(3); % 0.5;
    else
        half_mult = 0.707;
    end

    minus_0_9s_value = minus_1_1s_value *  inhibition_dist(4);
    minus_0_5s_value = minus_1_1s_value * half_mult;
    minus_0_1s_value = minus_1_1s_value * inhibition_dist(2);

    plus_0_9s_value = plus_1_1s_value * inhibition_dist(4);
    plus_0_5s_value = plus_1_1s_value * half_mult;
    plus_0_1s_value = plus_1_1s_value * inhibition_dist(2);

    % Print the parameter values
    [plus_1s_value plus_1_1s_value plus_0_5s_value minus_1s_value minus_1_1s_value minus_0_5s_value plus_3s_value];
    con_matrix(plus_1s_ind)    = plus_1s_value;con_matrix(plus_3s_ind)    = plus_3s_value;con_matrix(plus_1_1s_ind)  = plus_1_1s_value;con_matrix(plus_0_9s_ind)  = plus_0_9s_value;con_matrix(plus_0_5s_ind)  = plus_0_5s_value;con_matrix(plus_0_1s_ind)  = plus_0_1s_value;con_matrix(minus_1s_ind)   = minus_1s_value;con_matrix(minus_1_1s_ind) = minus_1_1s_value;con_matrix(minus_0_9s_ind) = minus_0_9s_value;con_matrix(minus_0_5s_ind) = minus_0_5s_value;con_matrix(minus_0_1s_ind) = minus_0_1s_value;

    % Add noise to synaptic weights
    % Alter synaptic weights of all synapses
    if strcmp(changed_synapse_type, 'All')
        % Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix)) .* con_matrix * (change_percent / 100);
        con_matrix = con_matrix + added_noise;
    end
    % Alter synaptic weights of one type of synapses on one side or randomly on both sides
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_left')
        con_matrix(P_EN_E_PG_output_synapses_indx_left) = con_matrix(P_EN_E_PG_output_synapses_indx_left) + con_matrix(P_EN_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_right')
        con_matrix(P_EN_E_PG_output_synapses_indx_right) = con_matrix(P_EN_E_PG_output_synapses_indx_right) + con_matrix(P_EN_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EN_E_PG_output_synapses_indx))) .* con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_alleq')
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_left')
        con_matrix(P_EG_E_PG_output_synapses_indx_left) = con_matrix(P_EG_E_PG_output_synapses_indx_left) + con_matrix(P_EG_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_right')
        con_matrix(P_EG_E_PG_output_synapses_indx_right) = con_matrix(P_EG_E_PG_output_synapses_indx_right) + con_matrix(P_EG_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EG_E_PG_output_synapses_indx))) .* con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_alleq')
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_left')
        con_matrix(E_PG_P_EN_output_synapses_indx_left) = con_matrix(E_PG_P_EN_output_synapses_indx_left) + con_matrix(E_PG_P_EN_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_right')
        con_matrix(E_PG_P_EN_output_synapses_indx_right) = con_matrix(E_PG_P_EN_output_synapses_indx_right) + con_matrix(E_PG_P_EN_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EN_output_synapses_indx))) .* con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_alleq')
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_left')
        con_matrix(E_PG_P_EG_output_synapses_indx_left) = con_matrix(E_PG_P_EG_output_synapses_indx_left) + con_matrix(E_PG_P_EG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_right')
        con_matrix(E_PG_P_EG_output_synapses_indx_right) = con_matrix(E_PG_P_EG_output_synapses_indx_right) + con_matrix(E_PG_P_EG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EG_output_synapses_indx))) .* con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_alleq')
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_left')
        con_matrix(E_PG_Pintr_output_synapses_indx_left) = con_matrix(E_PG_Pintr_output_synapses_indx_left) + con_matrix(E_PG_Pintr_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_right')
        con_matrix(E_PG_Pintr_output_synapses_indx_right) = con_matrix(E_PG_Pintr_output_synapses_indx_right) + con_matrix(E_PG_Pintr_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_Pintr_output_synapses_indx))) .* con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_alleq')
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
    end
    
    % Update the maxSpikeRate values with the ones read from the file
    if isKey(StimulusParamsMap, 'selectedColumnsBool')
        paramsMap = StimulusParamsMap('selectedColumnsBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'jumpBarBool')
        paramsMap = StimulusParamsMap('jumpBarBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'ImbalanceRightBool')
        paramsMap = StimulusParamsMap('ImbalanceRightBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesFixedBool')
        paramsMap = StimulusParamsMap('vonMisesFixedBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesTrajectoryBool')
        paramsMap = StimulusParamsMap('vonMisesTrajectoryBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    
    % Run the network simulation and get bump width at the end
    bump = runNetworkSimulation(con_matrix, display_plots, fixedThreshold, currentStrenth, jumpDistance, pulseLengthMulti, maxSpikeRate, inputList, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, SimulationParamsStruct);
    bumpWidth = bump.bumpWidth;
    % Print the bump width
    % bumpWidth * 180 / pi

    % bump.bumpLocPre
    % bump.bumpLocPost

    % bump.bumpWidthPre  * 180 / pi
    % bump.bumpWidthPost * 180 / pi

    % bump.spikeRate % Display the spike rates of the neurons
    
    out = bump;
    out.n_P_EN  = n_P_EN;
    out.n_P_EG  = n_P_EG;
    out.n_E_PG  = n_E_PG;
    out.n_Pintr = n_Pintr;
    out.con_matrix_filename = con_matrix_filename;
    out.con_matrix_parameter_set = parameter_set_locust;
    out.inhibition_distr_type = inhibition_distr_type;
    out.inhibition_width_sigma = inhibition_width_sigma;
    out.changed_synapse_type = changed_synapse_type;
    out.change_percent = change_percent;
    out.con_matrix = con_matrix;
    out.inputList = inputList;
    out.maxSpikeRate = maxSpikeRate;
    
end


%% Mixed Species connectivity With P-EG neurons
%% It uses the connectivity of Drosophila with inhibitory connections
%% replaced with those of the locust model.
function out=mixed_species_dros_with_noise(changed_synapse_type, change_percent, display_plots, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, jumpDistance, SimulationParamsStruct)
    fixedThreshold = false;
    jumpDistance = jumpDistance;
    currentStrenth = 0.0;
    pulseLengthMulti = 3;
    maxSpikeRate = 170; % 170 Hz is the original value used.

    % Mixed Species Drosophila corrected (2018) connectivity matrix including P-EG neurons
    con_matrix_filename = 'connectivity_matrix_drosophila_mine_case_5_9cols_localised_3_labels1.mat';
    load(con_matrix_filename);
    disp('Mixed Species Drosophila corrected (2018) connectivity matrix including P-EG neurons')
    
    if strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'E-PG')
        %           [       35      :       52           ];
        inputList = [n_P_EN+n_P_EG+1:n_P_EN+n_P_EG+n_E_PG]; % E-PG neurons indexes for Drosophila with 9 columns connected
    elseif strcmp(StimulusParamsMap('stimulatedNeuronGroup'), 'P-EN')
        %           [        1      :       16           ];
        inputList = [1:n_P_EN]; % P-EN neurons indexes for Drosophila with 9 columns connected
    end
    
    % Find matrix elements with these label values
    plus_1s_ind    = con_matrix == 1;plus_3s_ind    = con_matrix == 3;plus_1_1s_ind  = con_matrix == 1.1;
    plus_0_9s_ind = con_matrix == 0.9;
    plus_0_5s_ind  = con_matrix == 0.5;
    plus_0_1s_ind = con_matrix == 0.1;
    minus_1s_ind   = con_matrix == -1;minus_1_1s_ind = con_matrix == -1.1;minus_0_5s_ind  = con_matrix == -0.5;
    minus_0_9s_ind = con_matrix == -0.9;
    minus_0_1s_ind = con_matrix == -0.1;

    % Put ones on the P_EN output synaptic value indexes
    base_n_x = 0;       % Start index
    range_n_x = n_P_EN; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EN_E_PG_output_synapses_indx,P_EN_E_PG_output_synapses_indx_left,P_EN_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);
    
    base_n_x = n_P_EN;  % Start index
    range_n_x = n_P_EG; % How many indeces from start
    base_n_y = n_P_EN+n_P_EG; % Start column index
    range_n_y = n_E_PG;     % How many columnsfrom start
    [P_EG_E_PG_output_synapses_indx,P_EG_E_PG_output_synapses_indx_left,P_EG_E_PG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = 0;               % Start index
    range_n_y = n_P_EN;         % How many indeces from start
    [E_PG_P_EN_output_synapses_indx,E_PG_P_EN_output_synapses_indx_left,E_PG_P_EN_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG; % Start index
    range_n_x = n_E_PG;         % How many indeces from start
    base_n_y = n_P_EN;          % Start index
    range_n_y = n_P_EG;         % How many indeces from start
    [E_PG_P_EG_output_synapses_indx,E_PG_P_EG_output_synapses_indx_left,E_PG_P_EG_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    base_n_x = n_P_EN + n_P_EG;  % Start index
    range_n_x = n_E_PG;          % How many indeces from start
    base_n_y = n_P_EN + n_P_EG + n_E_PG; % Start index
    range_n_y = n_Pintr;                 % How many indeces from start
    [E_PG_Pintr_output_synapses_indx,E_PG_Pintr_output_synapses_indx_left,E_PG_Pintr_output_synapses_indx_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y);

    % <entry version> <fitness> <bumpLocPre> <bumpLocPost> <bumpWidthPre> <final bump width> <E-PG peak spike rate> <5 or 8 x weight parameters>
    inhibition_distr_type = 'uniform';
    inhibition_width_sigma = 0.0;
    parameter_set_drosophila = [5	0.0477273933933868	1.50000000000000	7	0.785398163397448	1.57079632679490	230	7.87200000000000	46.4385488157865	-53.6715878469669	-9.08917570780942	100.004200000000];
    
    if size(parameter_set_drosophila, 2) == 5
        [plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
    elseif size(parameter_set_drosophila, 2) == 6
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        plus_3s_value    = plus_1s_value * 3;
    elseif size(parameter_set_drosophila, 2) == 7
        [plus_1s_value,plus_1_1s_value,plus_0_5s_value,minus_1s_value,minus_1_1s_value,minus_0_5s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
    elseif size(parameter_set_drosophila, 2) == 9
        [version,fitness,bump_width,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        maxSpikeRate = peak_spike_rate_E_PG;
    elseif size(parameter_set_drosophila, 2) == 12
        [version,fitness,bump_loc_pre,bump_loc_post,bump_width_pre,bump_width_post,peak_spike_rate_E_PG,plus_1s_value,plus_1_1s_value,minus_1s_value,minus_1_1s_value,plus_3s_value] = feval(@(x) x{:}, num2cell(parameter_set_drosophila));
        bump_width = bump_width_post;
        maxSpikeRate = peak_spike_rate_E_PG;
    end

    % If the connectivity matrix contains -0.9 values then is sinusoidal
    % otherwise is the hardwired ratio
    if any(any(minus_0_9s_ind))
        half_mult = 0.5;
    else
        half_mult = 0.707;
    end

    minus_0_9s_value = minus_1_1s_value * 0.85355339;
    minus_0_5s_value = minus_1_1s_value * half_mult;
    minus_0_1s_value = minus_1_1s_value * 0.14644661;

    plus_0_9s_value = plus_1_1s_value * 0.85355339;
    plus_0_5s_value = plus_1_1s_value * half_mult;
    plus_0_1s_value = plus_1_1s_value * 0.14644661;

    % Print the parameter values
    [plus_1s_value plus_1_1s_value plus_0_5s_value minus_1s_value minus_1_1s_value minus_0_5s_value plus_3s_value];
    con_matrix(plus_1s_ind)    = plus_1s_value;con_matrix(plus_3s_ind)    = plus_3s_value;con_matrix(plus_1_1s_ind)  = plus_1_1s_value;con_matrix(plus_0_9s_ind)  = plus_0_9s_value;con_matrix(plus_0_5s_ind)  = plus_0_5s_value;con_matrix(plus_0_1s_ind)  = plus_0_1s_value;con_matrix(minus_1s_ind)   = minus_1s_value;con_matrix(minus_1_1s_ind) = minus_1_1s_value;con_matrix(minus_0_9s_ind) = minus_0_9s_value;con_matrix(minus_0_5s_ind) = minus_0_5s_value;con_matrix(minus_0_1s_ind) = minus_0_1s_value;

    % Add noise to synaptic weights
    % Alter synaptic weights of all synapses
    if strcmp(changed_synapse_type, 'All')
        % Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix)) .* con_matrix * (change_percent / 100);
        con_matrix = con_matrix + added_noise;
    end
    % Alter synaptic weights of one type of synapses on one side or randomly on both sides
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_left')
        con_matrix(P_EN_E_PG_output_synapses_indx_left) = con_matrix(P_EN_E_PG_output_synapses_indx_left) + con_matrix(P_EN_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_right')
        con_matrix(P_EN_E_PG_output_synapses_indx_right) = con_matrix(P_EN_E_PG_output_synapses_indx_right) + con_matrix(P_EN_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EN_E_PG_output_synapses_indx))) .* con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EN_to_E_PG_alleq')
        con_matrix(P_EN_E_PG_output_synapses_indx) = con_matrix(P_EN_E_PG_output_synapses_indx) + con_matrix(P_EN_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_left')
        con_matrix(P_EG_E_PG_output_synapses_indx_left) = con_matrix(P_EG_E_PG_output_synapses_indx_left) + con_matrix(P_EG_E_PG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_right')
        con_matrix(P_EG_E_PG_output_synapses_indx_right) = con_matrix(P_EG_E_PG_output_synapses_indx_right) + con_matrix(P_EG_E_PG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(P_EG_E_PG_output_synapses_indx))) .* con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'P_EG_to_E_PG_alleq')
        con_matrix(P_EG_E_PG_output_synapses_indx) = con_matrix(P_EG_E_PG_output_synapses_indx) + con_matrix(P_EG_E_PG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_left')
        con_matrix(E_PG_P_EN_output_synapses_indx_left) = con_matrix(E_PG_P_EN_output_synapses_indx_left) + con_matrix(E_PG_P_EN_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_right')
        con_matrix(E_PG_P_EN_output_synapses_indx_right) = con_matrix(E_PG_P_EN_output_synapses_indx_right) + con_matrix(E_PG_P_EN_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EN_output_synapses_indx))) .* con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EN_alleq')
        con_matrix(E_PG_P_EN_output_synapses_indx) = con_matrix(E_PG_P_EN_output_synapses_indx) + con_matrix(E_PG_P_EN_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_left')
        con_matrix(E_PG_P_EG_output_synapses_indx_left) = con_matrix(E_PG_P_EG_output_synapses_indx_left) + con_matrix(E_PG_P_EG_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_right')
        con_matrix(E_PG_P_EG_output_synapses_indx_right) = con_matrix(E_PG_P_EG_output_synapses_indx_right) + con_matrix(E_PG_P_EG_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_P_EG_output_synapses_indx))) .* con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_P_EG_alleq')
        con_matrix(E_PG_P_EG_output_synapses_indx) = con_matrix(E_PG_P_EG_output_synapses_indx) + con_matrix(E_PG_P_EG_output_synapses_indx) * (change_percent / 100);
    end
    
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_left')
        con_matrix(E_PG_Pintr_output_synapses_indx_left) = con_matrix(E_PG_Pintr_output_synapses_indx_left) + con_matrix(E_PG_Pintr_output_synapses_indx_left) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_right')
        con_matrix(E_PG_Pintr_output_synapses_indx_right) = con_matrix(E_PG_Pintr_output_synapses_indx_right) + con_matrix(E_PG_Pintr_output_synapses_indx_right) * (change_percent / 100);
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_all')
        %             Gaussian distributed mu=0 sigma=1 * value * percent change
        added_noise = randn(size(con_matrix(E_PG_Pintr_output_synapses_indx))) .* con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + added_noise;
    end
    if strcmp(changed_synapse_type, 'E_PG_to_Pintr_alleq')
        con_matrix(E_PG_Pintr_output_synapses_indx) = con_matrix(E_PG_Pintr_output_synapses_indx) + con_matrix(E_PG_Pintr_output_synapses_indx) * (change_percent / 100);
    end
    
    % Update the maxSpikeRate values with the ones read from the file
    if isKey(StimulusParamsMap, 'selectedColumnsBool')
        paramsMap = StimulusParamsMap('selectedColumnsBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'jumpBarBool')
        paramsMap = StimulusParamsMap('jumpBarBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'ImbalanceRightBool')
        paramsMap = StimulusParamsMap('ImbalanceRightBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesFixedBool')
        paramsMap = StimulusParamsMap('vonMisesFixedBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    if isKey(StimulusParamsMap, 'vonMisesTrajectoryBool')
        paramsMap = StimulusParamsMap('vonMisesTrajectoryBool');
        paramsMap('maxSpikeRate') = maxSpikeRate;
    end
    
    % Run the network simulation and get bump width at the end
    bump = runNetworkSimulation(con_matrix, display_plots, fixedThreshold, currentStrenth, jumpDistance, pulseLengthMulti, maxSpikeRate, inputList, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, SimulationParamsStruct);
    bumpWidth = bump.bumpWidth;
    % Print the bump width
    %bumpWidth * 180 / pi

    %bump.bumpLocPre
    %bump.bumpLocPost

    %bump.bumpWidthPre  * 180 / pi
    %bump.bumpWidthPost * 180 / pi

    %bump.spikeRate % Display the spike rates of the neurons

    out = bump;
    out.n_P_EN  = n_P_EN;
    out.n_P_EG  = n_P_EG;
    out.n_E_PG  = n_E_PG;
    out.n_Pintr = n_Pintr;
    out.con_matrix_filename = con_matrix_filename;
    out.con_matrix_parameter_set = parameter_set_drosophila;
    out.inhibition_distr_type = inhibition_distr_type;
    out.inhibition_width_sigma = inhibition_width_sigma;
    out.changed_synapse_type = changed_synapse_type;
    out.change_percent = change_percent;
    out.con_matrix = con_matrix;
    out.inputList = inputList;
    out.maxSpikeRate = maxSpikeRate;
    
end



%% Returns a matrix of the same dimensions as the input one with marked 
%% with logic 1's only the cells in the specified region that have values
%% other than 0. 
function [ret_matrix,ret_matrix_left,ret_matrix_right]=return_logic_marked_matrix_v2(con_matrix, base_n_x, range_n_x, base_n_y, range_n_y)
    ret_matrix = zeros(size(con_matrix));
    matrix_subset = con_matrix(base_n_y+1:base_n_y+range_n_y,base_n_x+1:base_n_x+range_n_x);
    matrix_subset(matrix_subset~=0) = true;
    ret_matrix(base_n_y+1:base_n_y+range_n_y,base_n_x+1:base_n_x+range_n_x) = matrix_subset;
    ret_matrix = logical(ret_matrix);
    
    ret_matrix_left = zeros(size(con_matrix));
    matrix_subset = con_matrix(base_n_y+1:base_n_y+range_n_y,base_n_x+1:floor(base_n_x+range_n_x/2));
    matrix_subset(matrix_subset~=0) = true;
    ret_matrix_left(base_n_y+1:base_n_y+range_n_y,base_n_x+1:floor(base_n_x+range_n_x/2)) = matrix_subset;
    ret_matrix_left = logical(ret_matrix_left);

    ret_matrix_right = zeros(size(con_matrix));
    matrix_subset = con_matrix(base_n_y+1:base_n_y+range_n_y,floor(base_n_x+range_n_x/2+1):base_n_x+range_n_x);
    matrix_subset(matrix_subset~=0) = true;
    ret_matrix_right(base_n_y+1:base_n_y+range_n_y,floor(base_n_x+range_n_x/2+1):base_n_x+range_n_x) = matrix_subset;
    ret_matrix_right = logical(ret_matrix_right);
end



% Plot the connectivity weight matrix
function plot_con_matrix(con_matrix, title_str, xlabel_str, ylabel_str, z_range_min, z_range_max)
    figure();
    surf(con_matrix);
    view(0, 90); % azimuth, elevation
    colorbar;     % show scale colorbar
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    title(title_str);
    % Set the range of values (colors)
    if ~exist('z_range_min', 'var')
        z_range_min = -40;
    end
    if ~exist('z_range_max', 'var')
        z_range_max = 40;
    end
    caxis([z_range_min z_range_max]);
    
    %% Set the zero value to correspond to black color in the colormap
    cmap = colormap; % Get the colormap values array
    cmap_portion = z_range_max / (z_range_max - z_range_min); % in [0, 1]
    indx = round(cmap_portion * size(cmap, 1)); % Index into the colormap corresponding to data value 0
    cmap(indx, :) = [0.5 0.5 0.5];
    cmap(indx-1, :) = [0.5 0.5 0.5];
    cmap(indx+1, :) = [0.5 0.5 0.5];
    colormap(cmap); % set the new array as the colormap
end

% Run the ring attractor simulation
function bump=runNetworkSimulation(con_matrix, display_plot_bool, fixedThreshold, currentStrenth, jumpDistance, pulseLengthMulti, maxSpikeRate, inputList, spiking_or_rate_based, StimulusParamsMap, initial_heading_index, SimulationParamsStruct)
    
    if isKey(StimulusParamsMap, 'simulationDuration')
        simulationTime = StimulusParamsMap('simulationDuration'); % in sec
    else
        simulationTime=60.0; % in sec
    end
    
    % stimulusSwitch=[0,0,0,0,0,0,0,0,0,0,0,0,1]; % Use the single pulse stimulus in specific neurons
    % New stimulsu specification method
    stimulusSwitch=zeros(1, 15);
    % Switch to new parameter passing and read values from 
    % the StimulusParamsMap map and override other arguments
    stimulusSwitch(14) = 1;
    % Do we apply stimulus to select neurons? 'selectedColumnsBool'
    if isKey(StimulusParamsMap, 'selectedColumnsBool')
        MP = StimulusParamsMap('selectedColumnsBool');
        stim_neuron_list = MP('stim_neuron_list');
        stimulusSwitch(13) = 1;
    end
    
    % Do we apply two pulses stimuli? 'jumpBarBool'
    if isKey(StimulusParamsMap, 'jumpBarBool')
        stimulusSwitch(10) = 1;
        stimulusSwitch(14) = 0;
    end
    
    % Do we apply uni-hemispheric stimulus? 'ImbalanceRightBool'
    if isKey(StimulusParamsMap, 'ImbalanceRightBool')
        stimulusSwitch(12) = 1;
    end
    
    % Do we apply von Mises stimulus to select neurons?
    if isKey(StimulusParamsMap, 'vonMisesFixedBool')
        stimulusSwitch(15) = 1;
    end
    
    % Do we apply a moving heading trajectory as a von Mises stimulus to neurons?
    if isKey(StimulusParamsMap, 'vonMisesTrajectoryBool')
        stimulusSwitch(16) = 1;
    end

    backgroundNoiseBool = 1;   % Background noise
    secondPulseTime = 0.25;    % Dummy for all cases but used by the jumpBarBool simulations
    initial_heading_index = initial_heading_index; % Dummy for all cases but used by the jumpBarBool simulations
    
    % If not specified use default value
    % It can be a list of neurons to stimulate on both hemispheres eg
    % stim_neuron_list = [4, 5];
    if ~exist('stim_neuron_list', 'var') || isempty(stim_neuron_list)
        stim_neuron_list = [floor(size(inputList,2)/2/2)]; % [4] for either 18 or 16 neurons
    end
    
    % Spiking or rate based neurons
    if ~exist('spiking_or_rate_based', 'var')
        spiking_or_rate_based = 1;
    else
        if strcmp(spiking_or_rate_based, 'spiking')
            spiking_or_rate_based = 1;
        elseif strcmp(spiking_or_rate_based, 'ratebased')
            spiking_or_rate_based = 3;
        elseif strcmp(spiking_or_rate_based, 'spiking_mem_gaussian_noise')
            spiking_or_rate_based = 4;
        end
    end
    
    % Run the ring attractor simulation
    out = PBexperimentModified(con_matrix', display_plot_bool, spiking_or_rate_based, currentStrenth, jumpDistance, pulseLengthMulti, simulationTime, stimulusSwitch, backgroundNoiseBool, inputList, maxSpikeRate, secondPulseTime, initial_heading_index, stim_neuron_list, StimulusParamsMap, SimulationParamsStruct);

    % Get the above threshold points with bluring
    blurredImage = out.blurredImage;

    % Get only the first 8 columns which must contain a bump
    firstColumns=blurredImage(:, 1:8);

    % Detect if there was a bump at the end of the simulation
    % Measure activity bump width as width at 50% of maximum
    lastActivity = firstColumns(end-2000:end-1000,:); % Get a segment in time
    lastActivity = mean(lastActivity, 1); % Get the mean of the segment
    bump2 = measureBump1(lastActivity, fixedThreshold);

    lastActivity = firstColumns(end-26000:end-25000,:); % Get a segment in time
    lastActivity = mean(lastActivity, 1); % Get the mean of the segment
    bump1 = measureBump1(lastActivity, fixedThreshold);
    
    % Set the result
    bump = bump2;
    bump.bumpWidthPre  = bump1.bumpWidth;
    bump.bumpLocPre    = bump1.bumpLoc;
    bump.bumpWidthPost = bump2.bumpWidth;
    bump.bumpLocPost   = bump2.bumpLoc;

    % Return the spike rates of neurons at the end of simulation
    bump.spikeRate = out.spikeRate;
    
    % Return also all other data
    bump.out = out;
    bump.simulationTime = simulationTime;    
    bump.stim_neuron_list = stim_neuron_list;
    
    bump3 = measureBump3(firstColumns, fixedThreshold, false);
    bump.r_ts = bump3.r_ts;
    bump.theta_ts = bump3.theta_ts;

end


