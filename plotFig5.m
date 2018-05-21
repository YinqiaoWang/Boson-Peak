dbstop if error
clearvars
close all

tic

quiverScale = 2;

m = 2;

pathHome = 'G:\YinqiaoWang\BosonPeak4\Data\BP4_1014';
folderList = {'180408a','180416a','180420a','180428a'};
nstepList = [10,7,6,1];

for jFolder = 1
    folderName = folderList{jFolder};
    
    path1 = [];
    folderListKt = [];
    load([pathHome,filesep,folderName,filesep,'path1.mat'],'path1','folderListKt');
    path1.cimg = [path1.data,filesep,'cimg'];
    
    path2 = [];
    load([path1.data,filesep,folderListKt{m},filesep,'path2.mat'],'path2');
    
    
    nstep = nstepList(jFolder);
    ksize = 1;

    eigenfrequency = [];
    eigenvector = [];
    fnt = [];
    load([path2.dynamicMatrix,filesep,sprintf('%06d',nstep),'.',sprintf('%02d',ksize),'.dynamicMatrix.mat'],...
        'eigenfrequency','eigenvector','diskID','fnt','FIS');
    eigenfrequency = eigenfrequency{1};
    eigenvector = eigenvector{1};
    fnt = fnt{1};
    diskID = diskID{1}(2:2:end)/2;
    
    centers = load([path1.centers,filesep,sprintf('%06d',nstep),'.center.txt']);

    idxFrequency = eigenfrequency > 720 & eigenfrequency < 1120;
    V2 = sum(eigenvector(:,idxFrequency).^2,2);
    V2 = sum(reshape(V2,2,[]),1);
    
    Ntop = 250;
    
    [V2_sort,idx_sort] = sort(V2,'descend');
    diskID_softspots = diskID(idx_sort(1:Ntop));
    
    %% Figure
    plotScale = 2;
    MarkerSize = plotScale*4;
    FontSize = plotScale*7;
    LineWidth = plotScale*1;

    cc = lines;
    
    P = {[0.08,0.6,0.38,0.38];
         [0.58,0.6,0.38,0.38];
         [0.08,0.1,0.38,0.38];
         [0.58,0.1,0.38,0.38]};
    
    
    figure
    set(gcf,'Units','centimeters','Position',[2,2,8,8]*plotScale)

    %% Coordination Number
    [binPdf1,edges] = histcounts(centers(diskID,5),'BinMethod','integers','Normalization','probability');
    [binPdf2,edges] = histcounts(centers(diskID_softspots,5),edges,'Normalization','probability');
    bins = (edges(1:end-1)+edges(2:end))/2;
    
    axes('Position',P{1})
    bar(bins,[binPdf1;binPdf2]')
    set(gca,'FontSize',FontSize,'xLim',[2,7],'XTick',[3,4,5,6],'LineWidth',LineWidth)
    xlabel('Coordination Number')
%     ylabel('Probability')
    legend({' All Particles',' Soft Spots'},'box','off')
    
    
    %% Local Inversion-Symmetry
    edges = linspace(0.4,1,6);
    [binPdf1] = histcounts(FIS,edges,'Normalization','probability');
    [binPdf2] = histcounts(FIS(idx_sort(1:Ntop)),edges,'Normalization','probability');
    bins = (edges(1:end-1)+edges(2:end))/2;
    
    axes('Position',P{2})
    bar(bins,[binPdf1;binPdf2]')
    set(gca,'FontSize',FontSize,'xLim',[0.4,1],'LineWidth',LineWidth)
    xlabel('$F_{IS}$','Interpret','latex')
%     ylabel('Probability')
%     legend({' All Particles',' Soft Spots'},'location','NorthWest')
    
    
    %% kn
    idx_soft_fnt = ismember(fnt(:,1),diskID_softspots);
    edges = linspace(0,4e3,7);
    [binPdf1] = histcounts(fnt(:,10),edges,'Normalization','probability');
    [binPdf2] = histcounts(fnt(idx_soft_fnt,10),edges,'Normalization','probability');
    bins = (edges(1:end-1)+edges(2:end))/2;
    
    axes('Position',P{3})
    bar(bins,[binPdf1;binPdf2]')
    set(gca,'FontSize',FontSize,'xLim',[0,4e3],'YLim',[0,0.4],'LineWidth',LineWidth)
    xlabel('$\langle k_n \rangle$','Interpret','latex')
%     ylabel('Probability')
%     legend({' All Particles',' Soft Spots'},'location','NorthEast')
    
    %% image
    img = imread([path2.Results,filesep,'LocalizedVibration.png']);
    axes('Position',P{4})
    imshow(img)
    
    saveas(gcf,[path2.Results,filesep,'Fig5.fig'])
    saveas(gcf,[path2.Results,filesep,'Fig5.jpg'])
    saveas(gcf,[path2.Results,filesep,'Fig5.pdf'])
    saveas(gcf,[path2.Results,filesep,'Fig5.eps'])
end
toc