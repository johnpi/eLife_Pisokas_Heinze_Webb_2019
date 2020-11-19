function out=flyLI2(params,connect,inputs,inputPSCs)

numNeurons=size(connect,1);
numInputs=length(inputs);

if numInputs ~=size(inputPSCs,2)
    error('input current count mismatch');
end

if isfield(params,'time');
    time=params.time; else time=4; end;                   % 0.5 (or 0.083s) seconds
if isfield(params,'dt');
    dt=params.dt; else dt=1e-4; end;                        % 1e-4 s
if isfield(params, 'N');
    N=params.N; else N=time/dt; end;
if isfield(params,'V0');
    %V0=params.V0; else V0=-52e-3; end;                      % -52e-3 (or -40e-3, -60e-3) V
    V0=params.V0; else V0=0; end;                      % -52e-3 (or -40e-3, -60e-3) V
if isfield(params,'Cmem');
    Cmem=params.Cmem; else Cmem=2e-8; end;                  %2e-7 (or 0.5e-7, 8e-7) F
if isfield(params,'Rmem');
    Rmem=params.Rmem; else Rmem=1e6; end;                   %1e6 Ohm

V=zeros(N,numNeurons)+V0;
V(1,:)=V0;
I=zeros(N,numNeurons);

% % Estimate spike rate with a low pass filter over stimulus spikes
% smoothing_interval = 100 * 10^-3; % [s] smothing window 100ms
% smoothing_samples_num = smoothing_interval / dt;
% stimulus1 = movmean(inputPSCs, smoothing_samples_num);
% % Scale up the stimulus values to the model. 
% stimulus1 = stimulus1 * 2100;


% Estimate spike rate with a low pass filter over stimulus spikes
smoothing_interval = 100 * 10^-3; % [s] smothing window 100ms
smoothing_samples_num = smoothing_interval / dt;
stimulus1 = movsum(inputPSCs, smoothing_samples_num) / smoothing_interval;
% Scale the stimulus spike rate values to the model. 
%stimulus1 = stimulus1 * 0.25;
stimulus1 = stimulus1 / max(max(stimulus1)) * 200; % Normalise to 200

% %stimulus1=zeros(N,16);
% stimulus1=zeros(N,18);
% pulseLength=0.25;  pulseLength=pulseLength/dt; %s
% pulseShift=pulseLength/2;
% experimentStartTime=0.5/dt;   %s
% %stimPhase=7;
% stimPhase=8;
% X=[0 0.5 1];
% Vq=[0 1 0];
% %for i=1:8
% for i=1:9
%     whereTemp=(pulseShift*(i-1)+1:pulseShift*(i-1)+pulseLength)+experimentStartTime;
%     XTemp=X*(whereTemp(end)-whereTemp(1))+whereTemp(1);
%     VTemp=interp1(XTemp,Vq,whereTemp);
%     %stimulus1(whereTemp,mod(i+stimPhase,8)+1)=VTemp;
%     stimulus1(whereTemp,mod(i+stimPhase,9)+1)=VTemp;
%     %stimulus1(whereTemp,mod(i+stimPhase,8)+1+8)=VTemp;
%     stimulus1(whereTemp,mod(i+stimPhase,9)+1+9)=VTemp;
% end
% disp('Leaky-integrator model is overwriting stimulus input configuration:');
% disp('Sweeping Bar v2 Stimulus used');
% 
% stim2Range=(2.25/dt):ceil(3.25/dt);
% stimulus1(stim2Range,[2 6 10 14])=1;
% disp('Competing Bar Stimulus used');

% Scale up stimulus1
%stimulus1 = stimulus1 * 100;

Iin=zeros(N,numNeurons);
dcCurrentScalingFactor=0.05; %determined empirically to give decent bump like parameters with default connectivity matrix
dcCurrentScalingFactor=0.06; % Working range 0.035 - 0.115
inCurrent=5e-9*dcCurrentScalingFactor;
Iin(:,inputs)=inCurrent*stimulus1;

currentNoise=0.15*2e-9; %optional current noise
%currentNoise=0.0; %optional current noise

for t=2:N
    for i=1:numNeurons
        %         size(connect(:,i))
        %         size((VtoI*(V(t-1,:)'-V0)))
        %         size(V)
        VtoI=V(t-1,:)'-V0;
        %VtoI=inCurrent*tanh(100*VtoI); % 100 determines the slope of the
        %sigmoid. inCurrent determines the amplitude.
        VtoI=inCurrent*tanh(100*VtoI);
        currentTemp=Iin(t-1,i) + sum(connect(:,i).*VtoI) +currentNoise*randn();
        dV=(1/Cmem)*( (V0-V(t-1,i))/Rmem + currentTemp );
        V(t,i)=V(t-1,i)+dV*dt;
    end
end

out.Iin=Iin;
out.V=V;
out.I=I;
out.connect=connect;
out.P=params;
out.P.time=time;
out.P.dt=dt;
out.P.V0=V0;
out.P.Cmem=Cmem;
out.P.Rmem=Rmem;
