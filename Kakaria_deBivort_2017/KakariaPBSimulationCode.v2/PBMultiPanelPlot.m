function figureOut=PBMultiPanelPlot(out)

    figure;
    Con=out.connect;
    blurredImage=out.blurredImage;
    P=out.P;

    molaspass=interp1([1 31 82 133 184 256],[0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
    blanding=interp1([1 128 129 256],[0 1 1; .248 .2559 .2559; .2559 .248 .248; 1 0 0],1:256);
    blanding=flipud(blanding);

    ax1=subplot(3,3,1);

    Con2=(Con-min(Con(:)))/(max(Con(:))-min(Con(:)));
    Con2=(Con2-mode(Con2(:)))*256+128;
    image(Con2)
    colormap(blanding)
    freezeColors

    ax2=subplot(3,3,4:9);

    imagesc(blurredImage');
    colormap(molaspass);
    hold on
    classBreaks=[8.5 16.5 24.5 32.5 41.5 50.5];
    for i=1:length(classBreaks)
        plot([0 P.N], classBreaks(i)*ones(1,2),'w-');
    end
    xlim([0 P.N]);
    colorbar
    
    ax3=subplot(3,3,1:3);
    hold on;
    tick=0;
    for i=1:size(out.inputPSCs,2);
        temp=out.inputPSCs(:,i);
%         temp=temp(2:end)-temp(1:end-1);
        temp=find(temp==1);
        scatter(temp,zeros(size(temp))-i,'k.');
    end
    
%     temp=out.Iin'>6*10^-8;
%     blurWidth=30*(1e-3/out.P.dt);
%     temp=imfilter(double(temp([33:40 43:50],:)),fspecial('gaussian',[1,blurWidth],blurWidth/5));
%     imagesc(temp);


    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
    
    linkaxes([ax1,ax2,ax3], 'x'); % link the x ases so zooming affects both
