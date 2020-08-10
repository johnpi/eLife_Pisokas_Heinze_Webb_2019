function out=flyLIF2(params,connect,inputs,inputPSCs,TRPCurrent)


numNeurons=size(connect,1);
numInputs=length(inputs);


% TRPCurrent=zeros(78,1);
% TRPCurrent(50:59)=pA*10^-9;

TRPCurrent=TRPCurrent*10^-9;

spikeThr=zeros(numNeurons,1)-45e-3;

if numInputs ~=size(inputPSCs,2)
    error('input current count mismatch');
end

if isfield(params,'time');
    time=params.time; else time=0.5; end;                   % 0.5 (or 0.083s) seconds
if isfield(params,'dt');
    dt=params.dt; else dt=1e-4; end;                        % 1e-4 s
if isfield(params, 'N');
    N=params.N; else N=time/dt; end;
if isfield(params,'V0');
    V0=params.V0; else V0=-52e-3; end;                      % -52e-3 (or -40e-3, -60e-3) V
if isfield(params,'spikeMin');
    spikeMin=params.spikeMin; else spikeMin=-72e-3; end;    %-72e-3 (or -90e-3) V
% if isfield(params,'spikeThr');
%     spikeThr=params.spikeThr; else spikeThr=-45e-3; end;    %-45e-3 (or -30e-3) V
if isfield(params,'Cmem');
    Cmem=params.Cmem; else Cmem=2e-8; end;                  %2e-8F
if isfield(params,'Rmem');
    Rmem=params.Rmem; else Rmem=1e6; end;                   %1e6 Ohm
if isfield(params,'apDuration');
    apDuration=params.apDuration; else apDuration=2; end;   %2 (or 0.5, 8) ms
if isfield(params,'pscRiseT');
    pscRiseT=params.pscRiseT; else pscRiseT=0.002; end;     %2 ms
if isfield(params,'pscFallT');
    pscFallT=params.pscFallT; else pscFallT=0.005; end;                 %5 ms
if isfield(params,'pscMag');
    pscMag=params.pscMag; else pscMag=5e-9; end;            %5e-9 (or 2.5e-9, 20e-9)A

if isfield(params,'pscWeights');
    pscWeights=params.pscWeights; else pscWeights=1; end;            %5e-9 (or 2.5e-9, 20e-9)A
if length(pscWeights)==1
    inputWeights=zeros(numInputs,1)+pscWeights;
end


spikeMax=20e-3;     %purely cosmetic parameter

apRiseTime=round((apDuration*10)/2);
apRise=normpdf(-1:1/apRiseTime:0);
apRise=(apRise-min(apRise))/(max(apRise)-min(apRise));
apRise=apRise*(spikeMax-spikeThr(1))+spikeThr(1);
apFallTime=round((apDuration*9)/2);
apFall=sin(pi()/2:pi()/apFallTime:3*pi()/2);
apFall=(apFall-min(apFall))/(max(apFall)-min(apFall));
apFall=apFall*(spikeMax-spikeMin)+spikeMin-.0001;
AP=[apRise apFall]';
AP(1)=[];
AP(end+1)=AP(end)+0.0001;
AP(end+1)=AP(end)+0.0001;
apLength=length(AP);

pscRise=sin(linspace(-pi()/2,pi()/2,pscRiseT/dt));
pscRise=(pscRise-min(pscRise))/(max(pscRise)-min(pscRise));
pscFall=2.^(-(0:(7*pscFallT/dt))*dt/pscFallT);
pscFall=(pscFall-min(pscFall))/(max(pscFall)-min(pscFall));
PSC=[pscRise pscFall]';
PSC=PSC*23.6305/sum(PSC); %normalize PSC magnitude to counteract changes in duration or time-constants
PSC=PSC*pscMag/max(PSC);

tBuff=length(PSC);

%prepare input current vectors
Iin=zeros(N+tBuff,numNeurons);
for i=1:numInputs
    if ~ iscell(inputPSCs)
        if min(inputPSCs(:))==0
        pscFallTimes=find(inputPSCs(:,i)==1);
        else
            pscFallTimes=inputPSCs(:,i);
        end
    else
        pscFallTimes=inputPSCs{i};
    end
    for j=1:length(pscFallTimes)
        PSCrange=pscFallTimes(j):(pscFallTimes(j)+tBuff-1);
        Iin(PSCrange,inputs(i))=Iin(PSCrange,inputs(i))+PSC*inputWeights(i);
    end
end

V=zeros(N+tBuff,numNeurons)+V0;
V(1,:)=V0;
I=zeros(N+tBuff,numNeurons);
APMask=zeros(N+tBuff,numNeurons);

for t=2:N
    for i=1:numNeurons
        if APMask(t,i)==0
            if V(t-1,i)>spikeThr(i)
                V(t-1,i)=spikeThr(i);
                V(t:t+apLength-1,i)=AP;
                APMask(t:t+apLength-1,i)=1;
                I(t:t+tBuff-1,i)=I(t:t+tBuff-1,i)+PSC;
            else
                currentTemp=Iin(t-1,i) + sum(connect(:,i).*I(t-1,:)') + TRPCurrent(i);
                %%currentTemp=Iin(t-1,i) + sum(connect(:,i).*I(t-1,:)') + TRPCurrent;
                dV=(1/Cmem)*( (V0-V(t-1,i))/Rmem + currentTemp );
                V(t,i)=V(t-1,i)+dV*dt;
            end
        end
    end
end

out.Iin=Iin;
out.V=V;
out.I=I;
out.inputPSCs=inputPSCs;
out.APMask=APMask;
out.connect=connect;
out.P=params;
out.P.time=time;
out.P.dt=dt;
out.P.V0=V0;
out.P.spikeMin=spikeMin;
out.P.spikeThr=spikeThr;
out.P.Cmem=Cmem;
out.P.Rmem=Rmem;
out.P.apDuration=apDuration;
out.P.pscRiseT=pscRiseT;
out.P.pscFallT=pscFallT;
