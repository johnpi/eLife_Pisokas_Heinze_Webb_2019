function out=PBexperiment(connect,visBool,lifMode,pA)

P.pscWeights=20;

if mod(sum(connect(:)),1)==0
    close all;
    numLabels=max(connect(:));
    %     signVector=(rand(numLabels,1)<0.5)*2-1;
    R1=(rand()>0.5)*2-1;
    R2=(rand()>0.5)*2-1;
    signVector=[1 1 1 1 1 -1 -1 -1 R1 R1 R1 R2 R2 R2 1 1 -1 -1];
    classMeans=(randn(numLabels,1)*P.pscWeights/3+P.pscWeights);
    classMeans=classMeans.*signVector'+abs(signVector')*1.5;
    classStds=rand(numLabels,1).*classMeans/100;
    
    conPack.connect=connect;
    conPack.means=classMeans;
    conPack.stds=classStds;
    Con=PBcon(conPack);
else
    Con=connect;
    classMeans='inherited';
    classStds='inherited';
end


P.time=4;              %4 or 2(s) for most figures
P.dt=1e-4;              %1e-4 by default
timeSteps=P.time/P.dt;  
P.N=timeSteps;

stimulus1=zeros(timeSteps,16);

pulseLength=0.25;  pulseLength=pulseLength/P.dt; %s, keep at 0.25
pulseShift=pulseLength/2;
sweep=1:pulseLength;
sweep=1-(sweep-pulseLength/2).^2;
sweep=1.5*(sweep-min(sweep))/(max(sweep)-min(sweep));
stimPhase=7;

baselineRate=5;    %Hz (5 in wowInstances)
sweepRate=120;       %Hz (120 in wowInstances

experimentStartTime=0.5/P.dt;   %s

SweepBarBool    =1; %=1 for multidimensional analysis
SweepBarBool2   =0;
BarsGoUpBool    =0;
DoubleSweepBool =0;
TwinBarBool     =0; %=1 for multidimensional analysis
ExtraInputBool  =0;
FatStripeBool   =0;
FlashingBarBool =0;
OneSpikeBool    =0;

P.inputDetails.SweepBarBool=SweepBarBool;
P.inputDetails.SweepBarBool2=SweepBarBool2;
P.inputDetails.BarsGoUpBool=BarsGoUpBool;
P.inputDetails.DoubleSweepBool=DoubleSweepBool;
P.inputDetails.TwinBarBool=TwinBarBool;
P.inputDetails.ExtraInputBool=ExtraInputBool;
P.inputDetails.FatStripeBool=FatStripeBool;
P.inputDetails.FlashingBarBool=FlashingBarBool;
P.inputDetails.OneSpikeBool=OneSpikeBool;

if SweepBarBool==1;
    for i=1:8
        whereTemp=(pulseShift*(i-1)+1:pulseShift*(i-1)+pulseLength)+experimentStartTime;
        stimulus1(whereTemp,mod(i+stimPhase,8)+1)=sweep;
        stimulus1(whereTemp,mod(i+stimPhase,8)+1+8)=sweep;
    end
end

if SweepBarBool2==1;
    X=[0 0.5 1];
    V=[0 1 0];
    for i=1:8
        whereTemp=(pulseShift*(i-1)+1:pulseShift*(i-1)+pulseLength)+experimentStartTime;
        XTemp=X*(whereTemp(end)-whereTemp(1))+whereTemp(1);
        VTemp=interp1(XTemp,V,whereTemp);
        stimulus1(whereTemp,mod(i+stimPhase,8)+1)=VTemp;
        stimulus1(whereTemp,mod(i+stimPhase,8)+1+8)=VTemp;
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
    disp('All possible inputs stimuluted');
    stimulus1=[stimulus1 zeros(timeSteps,4)];
    inputList=[33:41 42:5 59:62];
    error('check input dims');
else
    inputList=[33:40 43:50];
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

stimulus1=(stimulus1*sweepRate+baselineRate)*P.dt;
stimulus1=rand(size(stimulus1))<stimulus1;

switch lifMode
    case 0
        out=flyLIF(P,Con,inputList,stimulus1);
    case 1
        out=flyLIF2(P,Con,inputList,stimulus1,pA);
    case 2
        out=flyLI(P,Con);
        out.inputPSCs=[];
end

out.classMeans=classMeans;
out.classStds=classStds;

Spikes=out.V'>0;
blurWidth=120*(1e-3/out.P.dt);
blurredImage=imfilter(double(Spikes),fspecial('gaussian',[1,blurWidth],blurWidth/5));
blurredImage=blurredImage';
out.blurredImage=blurredImage;

if visBool==1
    PBMultiPanelPlot(out);
end


end