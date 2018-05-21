close all
clearvars
dbstop if error

tic

m = 2;

pathHome = 'G:\YinqiaoWang\BosonPeak4\Data\BP4_1014';
folderList = {'180408a','180416a','180420a','180428a'};
nstepList = [10,7,6,1];

for jFolder = 1
    folderName = folderList{jFolder};
    
    path1 = [];
    folderListKt = [];
    load([pathHome,filesep,folderName,filesep,'path1.mat'],'path1','folderListKt');
    
    path2 = [];
    load([path1.data,filesep,folderListKt{m},filesep,'path2.mat'],'path2');
    path2.modulus = [path2.data,filesep,'modulus'];
    
    nstep = nstepList(jFolder);
    ksize = 1;
    
    load([path2.dynamicStructureFactor,filesep,sprintf('%06d',nstep),'.',...
        sprintf('%02d',ksize),'.SCkwFitting2.mat'],...
        'wbins','kbins','CL','CT',...
        'modelfunCkLRZ','betaTCkLRZ','betaLCkLRZ',...
        'betaTCwLRZ','betaLCwLRZ')
    
    load([path2.dynamicStructureFactor,filesep,sprintf('%06d',nstep),'.',...
        sprintf('%02d',ksize),'.wmaxTL.mat'],'wmaxT','wmaxL','kbins')
    
    w0 = 2220.8;
    
    wBP = [720,1120]/w0;
    
    dw = 140;   % \delta \omega
    nbins = 50; % number of bins
    wbinEdges = dw*(0:nbins);
     
    wbin_counts = histcounts(wbins,wbinEdges);
    CL2 = cell2mat(cellfun(@(A) sum(A,1),mat2cell(CL',wbin_counts),'UniformOutput',false))*20/dw;
    CT2 = cell2mat(cellfun(@(A) sum(A,1),mat2cell(CT',wbin_counts),'UniformOutput',false))*20/dw;
    CL2 = CL2';
    CT2 = CT2';
    
    wbins2 = (wbinEdges(1:end-1)+wbinEdges(2:end))/2;
        
    plotScale = 2;
    MarkerSize = plotScale*3;
    FontSize = plotScale*7;
    LineWidth = plotScale*1;

    cc = lines;
    
    figure
    set(gcf,'Units','centimeters','Position',[2,2,8,8]*plotScale)

    %% Part1
    knnT = 3*(4:8);
    knnL = 1*(6:10);
    wlin = linspace(1,6000,10000);
    
    ax1 = axes('Position',[0.1 0.585 0.4 0.4]);
    hold on
    h = cell(numel(knnT),1);
    for jL = 1:numel(knnT)
        plot(wbins/w0,CT(knnT(jL),:),'s','LineStyle','none','MarkerSize',MarkerSize,...
            'Color',cc(jL,:));
    end
    legend({[repmat('$k = ',numel(knnT),1),num2str(kbins(knnT)','%01.2f'),repmat('$',numel(knnT),1)]},...
        'FontSize',FontSize,'Interpreter','latex','AutoUpdate','off','box','off',...
        'Position',[0.38,0.86,0.01,0.01])
    for jL = 1:numel(knnT)
        plot(wlin/w0,modelfunCkLRZ(betaTCkLRZ(knnT(jL),:),wlin),'-',...
            'Color',cc(jL,:),'LineWidth',LineWidth)
    end
    set(gca,'FontSize',FontSize,'LineWidth',LineWidth,...
        'XLim',[0.1,0.7],...
        'XMinorTick','on','box','on')
    ax1.XAxis.MinorTickValues = [0.3,0.5];
    
    xlabel('$\omega$','Interpret','latex','FontSize',FontSize)
    ylabel('$C_T(k,\omega)$','Interpret','latex','FontSize',FontSize)
   
    %% Part 2
    ax2 = axes('Position',[0.1 0.085 0.4 0.4]);
    hold on
    for jL = 1:numel(knnL)
        plot(wbins/w0,CL(knnL(jL),:),'o','LineStyle','none','MarkerSize',MarkerSize,...
            'Color',cc(jL,:),'MarkerFaceColor',cc(jL,:))
    end
    legend({[repmat('$k = ',numel(knnL),1),num2str(kbins(knnL)','%01.2f'),repmat('$',numel(knnL),1)]},...
        'FontSize',FontSize,'Interpreter','latex','AutoUpdate','off','box','off',...
        'Position',[0.38,0.36,0.01,0.01])
    for jL = 1:numel(knnL)
        plot(wlin/w0,modelfunCkLRZ(betaLCkLRZ(knnL(jL),:),wlin),'-',...
            'Color',cc(jL,:),'LineWidth',LineWidth)
    end
    set(gca,'FontSize',FontSize,'LineWidth',LineWidth,...
        'XLim',[0.1,0.7],...
        'XMinorTick','on','box','on')
    ax2.XAxis.MinorTickValues = [0.3,0.5];
    
    xlabel('$\omega$','Interpret','latex','FontSize',FontSize)
    ylabel('$C_L(k,\omega)$','Interpret','latex','FontSize',FontSize)
    
    %% Part 3
    [~,idx_kD] = min(abs(kbins-3.529)); 
    knn = idx_kD+8*(-3:4);
    
    axes('Position',[0.52 0.085 0.45 0.9],'Units','normalized')
    hold on
    for jL = 1:numel(knn)
        kn = knn(jL);
        Cmax = max(CL2(kn,:));
        plot(wbins2/w0,CL2(kn,:)/Cmax+jL-1,'o','LineStyle','none','MarkerSize',MarkerSize,...
            'Color',cc(jL,:),'MarkerFaceColor',cc(jL,:))

        plot(wbins2/w0,CT2(kn,:)/Cmax+jL-1,'s','LineStyle','none','MarkerSize',MarkerSize,...
            'Color',cc(jL,:))

%         text(0,jL-0.7,['$k = ',num2str(kbins(kn),'%01.2f'),'$'],...
%             'units','data','verticalalignment','top','FontSize',FontSize,'Interpret','latex',...
%             'Rotation',55)
    end
    set(gca,'FontSize',FontSize,'LineWidth',LineWidth,'XLim',[0,3],...
        'XTick',[0,1,2,3],'XMinorTick', 'on','yTickLabel',[],'YGrid','on','box','on')
    ax1 = gca;
    ax1.XAxis.MinorTickValues = (0.5:1:2.5);
    
    xlabel('$\omega$','Interpret','latex','FontSize',FontSize)
%     ylabel('$C(k,\omega)$','Interpret','latex','FontSize',FontSize)
    
    P = get(gca,'Position');
    axes('Position',P)
    hold on
    plot(wmaxL/w0,kbins,'k-','LineWidth',LineWidth)
    plot(wmaxT(kbins<3.5)/w0,kbins(kbins<3.5),'k-.','LineWidth',LineWidth)
    patch([wBP(1),wBP(2),wBP(2),wBP(1)],[0,0,10,10],'k','FaceAlpha',.3,'EdgeColor','none')
    hold off
    set(gca,'XLim',[0 3],'YLim',[2*kbins(knn(1))-kbins(knn(2)),kbins(knn(end))],'Visible','off')
    
    saveas(gcf,[path2.Results,filesep,'Fig2.fig'])
    saveas(gcf,[path2.Results,filesep,'Fig2.jpg'])
    saveas(gcf,[path2.Results,filesep,'Fig2.pdf'])
    saveas(gcf,[path2.Results,filesep,'Fig2.eps'])
end
toc