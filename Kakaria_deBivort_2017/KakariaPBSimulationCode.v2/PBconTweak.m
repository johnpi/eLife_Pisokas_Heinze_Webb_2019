function out=PBconTweak(con,baseCon,mode,params)

classes=unique(baseCon(baseCon~=0));
numTypes=length(classes);
conTemp=con;

switch(mode)
    case 'classNoise'
        whichClass=params(1);
        ci=params(2);
        classMean=mean(con(baseCon==whichClass));
        conTemp(baseCon==whichClass)=classMean+ci*classMean*randn(sum(baseCon(:)==whichClass),1);
    case 'noise'
        interClassNoise=params(1);
        intraClassNoise=params(2);
        for i=1:numTypes
            ci=mean(mean(con(baseCon==classes(i))));
            classOffset=ci*interClassNoise*randn();
            conTemp(baseCon==classes(i))=ci+classOffset+randn(sum(baseCon(:)==classes(i)),1)*intraClassNoise*ci;
        end
    case 'shift'
        if length(params)~=numTypes; error('length of shift vector does not match number of synapse types.'); end;
        for i=1:numTypes
            conTemp(baseCon==classes(i))=conTemp(baseCon==classes(i))*(1+params(i));
        end
end

out=conTemp;