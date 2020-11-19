%% secondPulseTime : 0-1 fraction of the simulation time at which the second input pulse is applied.
%% new_params_map  : A hash map of variable name to value, used for the new parameter passing style.
function out=PBexperimentModified(connect,visBool,lifMode,pA, jumpDistance, pulseLengthMulti, simulationTime, stimulusSwitch, backgroundNoiseBool, inputList, maxSpikeRate, secondPulseTime, initial_heading_index, stim_neuron_list, new_params_map, SimulationParamsStruct)

P.pscWeights=20;

Con=connect;
% If given a single value, typically 0, create a vector of 0.
if length(pA) == 1
    pA = pA * ones(1, size(Con,1));
end

classMeans='inherited';
classStds='inherited';

% Optional specification of the simulation time
if exist('simulationTime', 'var')
    P.time = simulationTime;
else
    P.time=4;              %4 or 2(s) for most figures
end
P.dt=1e-4;              %1e-4 by default
timeSteps=P.time/P.dt;  
timeSteps=floor(timeSteps); % added to allow simulationTime to result to integer
P.N=timeSteps;

% Optional variable to change the default duration of stimulus
if ~exist('pulseLengthMulti', 'var')
    pulseLengthMulti = 1;
end

pulseLength=0.25*pulseLengthMulti;  pulseLength=pulseLength/P.dt; %s, keep at 0.25
pulseShift=pulseLength/2;
sweep=1:pulseLength;
sweep=1-(sweep-pulseLength/2).^2;
sweep=1.5*(sweep-min(sweep))/(max(sweep)-min(sweep));
stimPhase=7;

if exist('inputList', 'var')
    cols = size(inputList, 2);
else
    cols = 16; % The default
end    
halfCols = cols / 2;      % 8 is the default value. For 18 columns will be 9
stimPhase = halfCols - 1; % 7 is the default value. For 18 columns will be 8

stimulus1=zeros(timeSteps,cols);


baselineRate=5;    %Hz (5 in wowInstances)

if ~exist('maxSpikeRate', 'var')
    sweepRate=120;       %Hz (120 in wowInstances
else
    sweepRate=maxSpikeRate; % Maximum firing rate in Hz
end

experimentStartTime=0.5/P.dt;   %s

if exist('SimulationParamsStruct', 'var')
    if isfield(SimulationParamsStruct, 'Cmem')
        P.Cmem = SimulationParamsStruct.Cmem; % Nominal membrane capacitance
    end
    if isfield(SimulationParamsStruct, 'Rmem')
        P.Rmem = SimulationParamsStruct.Rmem;  % Nominal membrane resistence
    end
    if isfield(SimulationParamsStruct, 'mem_noise')
        P.mem_noise = SimulationParamsStruct.mem_noise; % Type of membrane parameters noise
    end
    if isfield(SimulationParamsStruct, 'Csigma')
        P.Csigma = SimulationParamsStruct.Csigma;     % Variance if the normally distributed noise to membrane parameters
    end
    if isfield(SimulationParamsStruct, 'Rsigma')
        P.Rsigma = SimulationParamsStruct.Rsigma;     % Variance if the normally distributed noise to membrane parameters
    end
end

% If stimulusSwitch vector of Bool values was given use it
if exist('stimulusSwitch', 'var')
    vec_length = 16;
    if size(stimulusSwitch, 2) < vec_length
        stimulusSwitch(length(stimulusSwitch)+1:vec_length) = 0;
    end
    [SweepBarBool,SweepBarBool2,BarsGoUpBool,DoubleSweepBool,TwinBarBool,ExtraInputBool,FatStripeBool,FlashingBarBool,OneSpikeBool,JumpBarBool,initBarBool,ImbalanceRightBool,selectedColumnsBool,useNewInputMapValues,vonMisesFixedBool,vonMisesTrajectoryBool] = feval(@(x) x{:}, num2cell(stimulusSwitch));
else % else use in file settings
    SweepBarBool    =1; %=1 for multidimensional analysis
    SweepBarBool2   =0;
    BarsGoUpBool    =0;
    DoubleSweepBool =0;
    TwinBarBool     =0; %=1 for multidimensional analysis
    ExtraInputBool  =0;
    FatStripeBool   =0;
    FlashingBarBool =0;
    OneSpikeBool    =0;
    JumpBarBool     =0; % My jump response stimulus
    initBarBool     =0; % My initialisation with a pulse
    ImbalanceRightBool = 0; % My application of stimulus on the right side
    % My option for providing equal stimulus to a list of neurons on both 
    % hemispheres specified in variable stim_neuron_list
    selectedColumnsBool = 0; 
    useNewInputMapValues = 0; % Switch to using variable values provided in the map new_params_map. That was done because adding new arguments for new functionality ended up being messy.
    vonMisesFixedBool = 0; % Apply a von Misses dstributed stimulus (circular Gaussian)
    vonMisesTrajectoryBool = 0; % Apply a moving von Misses dstributed stimulus (circular Gaussian)
end

P.inputDetails.SweepBarBool=SweepBarBool;
P.inputDetails.SweepBarBool2=SweepBarBool2;
P.inputDetails.BarsGoUpBool=BarsGoUpBool;
P.inputDetails.DoubleSweepBool=DoubleSweepBool;
P.inputDetails.TwinBarBool=TwinBarBool;
P.inputDetails.ExtraInputBool=ExtraInputBool;
P.inputDetails.FatStripeBool=FatStripeBool;
P.inputDetails.FlashingBarBool=FlashingBarBool;
P.inputDetails.OneSpikeBool=OneSpikeBool;
P.inputDetails.JumpBarBool=JumpBarBool;
P.inputDetails.initBarBool = initBarBool;
P.inputDetails.ImbalanceRightBool = ImbalanceRightBool;
P.inputDetails.selectedColumnsBool = selectedColumnsBool;
P.inputDetails.useNewInputMapValues = useNewInputMapValues; % Switch to using variable values provided in the map new_params_map. That was done because adding new arguments for new functionality ended up being messy.
P.inputDetails.vonMisesFixedBool = vonMisesFixedBool;
P.inputDetails.vonMisesTrajectoryBool = vonMisesTrajectoryBool;

if ~exist('backgroundNoiseBool', 'var')
    backgroundNoiseBool = 1;
end
P.inputDetails.backgroundNoiseBool = backgroundNoiseBool;

if ~exist('secondPulseTime', 'var')
    secondPulseTime = 0.25;
end
if secondPulseTime > 1
    secondPulseTime = 0.25;
end

if ~exist('initial_heading_index', 'var')
    initial_heading_index = 2;
end

if selectedColumnsBool==1 && ~exist('stim_neuron_list', 'var')
    stim_neuron_list = [halfCols/2];
end

if selectedColumnsBool==1
    pulseStartTime = 0;
    if useNewInputMapValues==1 && isKey(new_params_map, 'selectedColumnsBool')
        % Override parameter values with those from the map
        map = new_params_map('selectedColumnsBool');
        pulseStartTime = map('pulseStartTime');
        pulseLength = map('pulseDuration');
        stim_neuron_list = map('stim_neuron_list');
        maxSpikeRateLocal = map('maxSpikeRate');
        maxSpikeRateFactorLocal = map('maxSpikeRateFactor');
        pulseStartTime=pulseStartTime/P.dt;
        pulseLength=pulseLength/P.dt;
        sweep=1:pulseLength;
        sweep=1-(sweep-pulseLength/2).^2;
        sweep=1.5*(sweep-min(sweep))/(max(sweep)-min(sweep));
        sweep = sweep * maxSpikeRateLocal * maxSpikeRateFactorLocal;
    end
    whereTemp=(pulseStartTime+1:pulseStartTime+pulseLength);
    % Assign the same sweep vector to all listed rows/neurons on both sides
    for col_indx = stim_neuron_list
        stimulus1(whereTemp,col_indx)=sweep;  % Left side
    end
    for col_indx = halfCols+stim_neuron_list
        stimulus1(whereTemp,col_indx)=sweep; % Right side
    end
end

if initBarBool==1
    % Initialise the ring attactor with input activity bump applied at
    % column 8/2 or floor(9/2)
    i=floor(halfCols/2);
        whereTemp=(1:pulseLength);
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1)=sweep;
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1+halfCols)=sweep;
end

if ImbalanceRightBool==1
    if useNewInputMapValues==1 && isKey(new_params_map, 'ImbalanceRightBool')
        % Override parameter values with those from the map
        map = new_params_map('ImbalanceRightBool');
        pulseStartTime = map('pulseStartTime');
        pulseLength = map('pulseDuration');
        stim_neuron_list = map('stim_neuron_list');
        maxSpikeRateLocal = map('maxSpikeRate');
        maxSpikeRateFactorLocal = map('maxSpikeRateFactor');
        pulseLength=pulseLength/P.dt;
        pulseStartTime=pulseStartTime/P.dt;
        % It was
        %sweep=1:pulseLength;
        %sweep=1-(sweep-pulseLength/2).^2;
        %sweep=1.5*(sweep-min(sweep))/(max(sweep)-min(sweep));
        % Make stimulus constant for the whole period
        sweep = ones(1,length(1:pulseLength));
        sweep = sweep * maxSpikeRateLocal * maxSpikeRateFactorLocal;
        % Apply one sided stimulation in all columns but one side
        whereTemp=(pulseStartTime+1:pulseStartTime+pulseLength);
        % Assign the same sweep vector to all rows on one side
        for col_indx = [1:halfCols]
            stimulus1(whereTemp,col_indx)=sweep;  % Left side
        end
        for col_indx = [halfCols+1:2*halfCols]
            %stimulus1(whereTemp,col_indx)=sweep; % Right side
        end
    else
        % Apply stronger input to the P-ENs on the right side of the PB
        secondPulseTimeStep = size(stimulus1, 1) * secondPulseTime;
        % Apply one sided stimulation in all columns but one side
        whereTemp=(secondPulseTimeStep+1:secondPulseTimeStep+pulseLength)+experimentStartTime;
        % Assign the same sweep vector to all rows on one side
        for col_indx = [1:halfCols]
            stimulus1(whereTemp,col_indx)=sweep;  % Left side
        end
        for col_indx = [halfCols+1:2*halfCols]
            %stimulus1(whereTemp,col_indx)=sweep; % Right side
        end
    end
end

if JumpBarBool==1
    % Moment of second input stimulus
    firstPulseTimeStep = 0;
    secondPulseTimeStep = size(stimulus1, 1) * secondPulseTime;
    % If angular jump distance is not given use default value
    if ~exist('jumpDistance', 'var')
        jumpDistance = 180; % Angular jump in degrees
    end
    %fprintf('jumpDistance = %f', jumpDistance);
    % Convert jump distance from degrees to column index=[1,8]
    initial_heading_index=initial_heading_index; % was =2;
    % How many neurons the stimulus must move
    add_heading_index = floor(mod(jumpDistance / 360 * halfCols, halfCols));
    % Add it to the initial position and convert it to wrap around the 8
    % neurons of the ring attractor
    second_heading_index=mod(initial_heading_index+add_heading_index+halfCols-1, halfCols) + 1;
    i=initial_heading_index; % Initial input activity bump applied at column 1
        whereTemp=(firstPulseTimeStep+1:firstPulseTimeStep+pulseLength) + 0; % was +experimentStartTime;
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1)=sweep;
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1+halfCols)=sweep;
    i=second_heading_index; % Second input activity bump at column heading_index
        whereTemp=(secondPulseTimeStep+1:secondPulseTimeStep+pulseLength)+experimentStartTime;
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1)=sweep;
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1+halfCols)=sweep;
end


if vonMisesFixedBool==1
    map = new_params_map('vonMisesFixedBool');
    pulseStartTime = map('pulseStartTime');
    pulseLength = map('pulseDuration');
    stim_neuron_list = map('stim_neuron_list');
    maxSpikeRateLocal = map('maxSpikeRate');
    maxSpikeRateFactorLocal = map('maxSpikeRateFactor');
    kappa = map('kappa');
    
    stim_neuron_center_deg = (stim_neuron_list(1) - 1) / halfCols * 360;
    %kappa = pi*7/6; % This results to Gaussian with 91deg FWHM
    sample_at_deg = [0:360/8:359];
    sample_at_rads = circ_ang2rad(sample_at_deg);
    [p alpha] = circ_vmpdf(sample_at_rads, circ_ang2rad(stim_neuron_center_deg), kappa);
    % Normalise p to the range of 0 to 1 because it varies with kappa
    p = (p - min(p)) * 1 / max(p - min(p));
    
    % In the drosophila with 9th columns copy over the value of the 1st
    % column
    if halfCols == 9
        p = [p; p(1)];
    end
    pulseStartTime=pulseStartTime/P.dt;
    pulseLength=pulseLength/P.dt;
    sweep=1:pulseLength;
    sweep=1-(sweep-pulseLength/2).^2;
    sweep=1.5*(sweep-min(sweep))/(max(sweep)-min(sweep));
    sweep = sweep * maxSpikeRateLocal * maxSpikeRateFactorLocal;

    whereTemp=(pulseStartTime+1:pulseStartTime+pulseLength);
    % Assign the same scaled sweep vector to all listed rows/neurons on both sides
    for col_indx = 1:halfCols
        stimulus1(whereTemp,col_indx)         =sweep*p(col_indx); % Left side
        stimulus1(whereTemp,col_indx+halfCols)=sweep*p(col_indx); % Right side
    end 
end

if vonMisesTrajectoryBool==1
    map = new_params_map('vonMisesTrajectoryBool');
    pulseStartTime = map('pulseStartTime');
    pulseLength = map('pulseDuration');
    %stim_neuron_list = map('stim_neuron_list'); % Nope
    stim_neuron_center_idx_ts = map('stim_neuron_list');
    maxSpikeRateLocal = map('maxSpikeRate');
    maxSpikeRateFactorLocal = map('maxSpikeRateFactor');
    kappa = map('kappa');
    
    stim_neuron_center_idx_ts_deg = (stim_neuron_center_idx_ts - 1) / halfCols * 360;
    %kappa = pi*7/6; % This results to Gaussian with 91deg FWHM
    sample_at_deg = [0:360/8:359];
    sample_at_rads = circ_ang2rad(sample_at_deg);
    [p alpha] = circ_vmpdf(sample_at_rads, circ_ang2rad(stim_neuron_center_idx_ts_deg), kappa);
    % Normalise p to the range of 0 to 1 because it varies with kappa
    p = (p - min(min(p))) * 1 / max(max(p - min(min(p))));
    
    % In the drosophila with 9th columns copy over the value of the 1st
    % column
    if halfCols == 9
        p = [p; p(1, :)];
    end
    pulseStartTime=pulseStartTime/P.dt;
    pulseLength=pulseLength/P.dt;
    sweep=ones(1, round(pulseLength));
    %sweep=1-(sweep-pulseLength/2).^2;
    %sweep=1.5*(sweep-min(sweep))/(max(sweep)-min(sweep));
    sweep = sweep * maxSpikeRateLocal * maxSpikeRateFactorLocal;

    whereTemp=(pulseStartTime+1:pulseStartTime+pulseLength);
    % Assign the same scaled sweep vector to all listed rows/neurons on both sides
    for col_indx = 1:halfCols
        stimulus1(whereTemp,col_indx)         =sweep.*p(col_indx,:); % Left side
        stimulus1(whereTemp,col_indx+halfCols)=sweep.*p(col_indx,:); % Right side
    end 
end



if SweepBarBool==1;
    for i=1:halfCols
        whereTemp=(pulseShift*(i-1)+1:pulseShift*(i-1)+pulseLength)+experimentStartTime;
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1)=sweep;
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1+halfCols)=sweep;
    end
end

if SweepBarBool2==1;
    X=[0 0.5 1];
    V=[0 1 0];
    for i=1:halfCols
        whereTemp=(pulseShift*(i-1)+1:pulseShift*(i-1)+pulseLength)+experimentStartTime;
        XTemp=X*(whereTemp(end)-whereTemp(1))+whereTemp(1);
        VTemp=interp1(XTemp,V,whereTemp);
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1)=VTemp;
        stimulus1(whereTemp,mod(i+stimPhase,halfCols)+1+halfCols)=VTemp;
    end
end


if BarsGoUpBool==1
    stimulus1=fliplr(stimulus1);
end

if DoubleSweepBool==1
    s1A=stimulus1(1:(1/P.dt),:);
    stimulus1((1/P.dt+1):(2/P.dt),:)=s1A;
end

stim2Range=(2.25/P.dt):ceil(3/P.dt);

if FatStripeBool==1
    stimulus1(stim2Range,[2:5 10:13])=1;
end

if TwinBarBool==1
    stimulus1(stim2Range,[2 6 10 14])=1;
end

if FlashingBarBool==1
    flashLength=100;
    flash=[zeros(flashLength,1) ones(flashLength,1)];
    flashes=[];
    stimRange=stim2Range;
    while length(flashes)<length(stimRange)
        flashes=[flashes flash];
    end
    flashes=flashes(1:length(stimRange)).*stimRange;
    flashes(flashes==0)=[];
    stimulus1(flashes,[4 12])=1;
end

if ExtraInputBool==1
    disp('All possible inputs stimulated');
    stimulus1=[stimulus1 zeros(timeSteps,4)];
    inputList=[33:41 42:5 59:62];
    error('check input dims');
else
    % inputList=[33:40 43:50];  % It was like this
    % Replaced with this
    if ~exist('inputList', 'var')
        inputList=[33:48];        % But my matrixes do not have the gap so this is now
    end
end

if OneSpikeBool==1
    stimulus1(1500,7)=inf;
    baselineRate=0;
end

P.inputDetails.stimPhase=stimPhase;
P.inputDetails.baselineRate=baselineRate;
P.inputDetails.sweepRate=sweepRate;
P.inputDetails.experimentStartTime=experimentStartTime;
P.inputDetails.pulseLength=pulseLength;
P.inputDetails.pulseShift=pulseShift;
P.inputDetails.inputList=inputList;

stimulus1_mask = stimulus1 > 0; % Store this in order to mask the spikes signal later
if useNewInputMapValues==0
    % Scale the sweep to the required peak spike rate and add base spike rate
    stimulus1=(stimulus1*sweepRate+baselineRate)*P.dt;
else
    % We have already scaled the sweep to have the required peak spike rate
    stimulus1=(stimulus1+baselineRate)*P.dt;
end
stimulus1=rand(size(stimulus1))<stimulus1; % Convert to spikes

% If this variable is false we should remove the background random spiking
if ~backgroundNoiseBool
    stimulus1 = and(stimulus1, stimulus1_mask);
end

switch lifMode
    case 0
        out=flyLIF(P,Con,inputList,stimulus1);
    case 1
        out=flyLIF2(P,Con,inputList,stimulus1,pA);
    case 2
        out=flyLI(P,Con);
        out.inputPSCs=[];
    case 3 % My implementation of rate based simulation
        out=flyLI2(P,Con,inputList,stimulus1);
        out.inputPSCs=[];        
    case 4 % My extension of spiking simulation with ability to modify individual neuron membrane properties
        % P can optionally have the properties P.Cmem and P.Rmem be, instead of scalar values, 
        % vectors of length equal to the number of neurons (side of connectivity matrix)
        if ~isfield(P, 'Cmem')
            P.Cmem = 2e-8; % Nominal membrane capacitance
        end
        if ~isfield(P, 'Rmem')
            P.Rmem = 1e6;  % Nominal membrane resistence
        end
        if ~isfield(P, 'Csigma')
            P.Csigma = 0;     % Variance of the normally distributed noise to membrane parameters
        end
        if ~isfield(P, 'Rsigma')
            P.Rsigma = 0;     % Variance of the normally distributed noise to membrane parameters
        end
        % P.mem_noise = 'gaussian';
        % Generate the vectors with the membrane properties of the neurons
        if isfield(P, 'mem_noise')
            if strcmp(P.mem_noise, 'gaussian')
                numNeurons=size(Con,1);
                P.Cmem = ones(1, numNeurons) * P.Cmem;
                P.Cmem = P.Cmem + normrnd(0, P.Csigma, size(P.Cmem)) .* P.Cmem;
                P.Cmem(P.Cmem<0) = 0; % Do not allow negative values
                P.Rmem = ones(1, numNeurons) * P.Rmem;
                P.Rmem = P.Rmem + normrnd(0, P.Rsigma, size(P.Rmem)) .* P.Rmem;
                P.Rmem(P.Rmem<0) = 0; % Do not allow negative values
            end
        end
        out=flyLIF2mod(P,Con,inputList,stimulus1,pA);
end

out.classMeans=classMeans;
out.classStds=classStds;

Spikes=out.V'>0;
blurWidth=120*(1e-3/out.P.dt);
blurredImage=imfilter(double(Spikes),fspecial('gaussian',[1,blurWidth],blurWidth/5));
blurredImage=blurredImage';
out.blurredImage=blurredImage;

% If using the rate based model we do not have spikes so it is as if there
% are constantly spikes resulting to high smoothed value and saturation
% for this reason use an alternative method. 
if lifMode == 3
    Vrectified = out.V;
    Vrectified(out.V<0) = 0;
    out.blurredImage = Vrectified;
    %out.blurredImage = out.blurredImage / max(max(out.blurredImage)); % normalise to max 1
    %out.blurredImage(out.blurredImage>0.7) = 0.7; % crop at max 0.7
end

%% My addition of measuring spike rates at the end of the simulation
out.spikeRate = countSpikeRateAtEnd(Spikes, P.dt);
%% And storing the spike sequence
out.spikeSequence = markSpikes(Spikes);

if visBool==1
    PBMultiPanelPlot(out);
end


end

%% My addition for calculating the spiking rate at the end of the simulation
function spikesSequence=markSpikes(SpikesArray)
    %% SpikesArray is assumed to have one row per neuron
    %% and one column for each time step.
    
    % Replace sequences of logical 1's representing single 
    % action potentials with a single 1's
    C = [];
    for i=[1:size(SpikesArray, 1)]
        A = SpikesArray(i, :);
        B=A;
        B(diff([0 A])==0)=0;
        C=[C; B];
    end
    
    spikesSequence = C;
end

%% My addition for calculating the spiking rate at the end of the simulation
function spikeRate=countSpikeRateAtEnd(SpikesArray, dt)
    %% SpikesArray is assumed to have one row per neuron
    %% and one column for each time step.
    
    % Replace sequences of logical 1's representing single 
    % action potentials with a single 1's
    C = markSpikes(SpikesArray);
    
    % Get a segment from the end of the time series (last 1000 time samples)
    c=C(:,end-999:end);
    s=sum(c, 2); % Count the 1's per row
    spikeRate = s ./ (size(c, 2) * dt); % Calculate 1's per sec
end

