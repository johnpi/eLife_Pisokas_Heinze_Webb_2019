function out=PBModulusMean(data)
nFrames=size(data,1);
nDims=size(data,2);
out=zeros(nFrames,1);

for i=1:nFrames
   dataTemp=[data(i,:);(1:nDims)];
   dataTemp=repmat(dataTemp,1,5);
   which=find(dataTemp(1,:)==max(dataTemp(1,:)));
   which=which(ceil(length(which)/2));

   dataTemp=dataTemp(:,(which-ceil(nDims/2)):(which-ceil(nDims/2)+nDims-1));
   centroid=sum((dataTemp(1,:).*(1:nDims)))/sum(dataTemp(1,:));
   out(i)=centroid+dataTemp(2,1)-1;
end

out=mod(out,nDims);
out(out<0.5)=out(out<0.5)+8;
