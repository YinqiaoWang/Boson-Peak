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
    
    load([path2.modulus,filesep,sprintf('%06d',nstep),'.',sprintf('%02d',ksize),'.modulus.mat'],'vL','vT');
    
    w0 = 2220.8;
    wBP = [720,1120]/w0;
    wbins = wbins/w0;
    
    vL0 = vL*10000/(51+70)/w0;    % m/s -> D/s
    vT0 = vT*10000/(51+70)/w0;    % m/s -> D/s
    
    plotScale = 2;
    MarkerSize = plotScale*4;
    FontSize = plotScale*7;
    LineWidth = plotScale*1;

    cc = lines;
    
    figure
    set(gcf,'Units','centimeters','Position',[2,2,16,6]*plotScale)
    
    %% Part1
    knnT = 3*(4:9);
    knnL = 1*(6:11);
    wlin = linspace(1,6000,10000);
    
    axes('Position',[0.07 0.15 0.25 0.8])
    hold on
    h = cell(numel(knnT),1);
    for jL = 1:numel(knnT)
        plot(wbins,CT(knnT(jL),:),'s','LineStyle','none','MarkerSize',MarkerSize,...
            'Color',cc(jL,:));
    end
    legend({[repmat('$k = ',numel(knnT),1),num2str(kbins(knnT)','%01.2f'),repmat('$',numel(knnT),1)]},...
        'FontSize',FontSize,'Interpreter','latex','AutoUpdate','off','box','off')
    for jL = 1:numel(knnT)
        plot(wlin/w0,modelfunCkLRZ(betaTCkLRZ(knnT(jL),:),wlin),'-',...
            'Color',cc(jL,:),'LineWidth',LineWidth)
    end
    set(gca,'FontSize',FontSize,'LineWidth',LineWidth,...
        'XLim',[0.1,0.7],...
        'XMinorTick','on','box','on')
    
    xlabel('$\omega$','Interpret','latex','FontSize',FontSize)
    ylabel('$C_T(k,\omega)$','Interpret','latex','FontSize',FontSize)
   
    %% Part 2
    axes('Position',[0.39 0.15 0.25 0.8])
    hold on
    for jL = 1:numel(knnL)
        plot(wbins,CL(knnL(jL),:),'o','LineStyle','none','MarkerSize',MarkerSize,...
            'Color',cc(jL,:),'MarkerFaceColor',cc(jL,:))
    end
    legend({[repmat('$k = ',numel(knnL),1),num2str(kbins(knnL)','%01.2f'),repmat('$',numel(knnL),1)]},...
        'FontSize',FontSize,'Interpreter','latex','AutoUpdate','off','box','off')
    for jL = 1:numel(knnL)
        plot(wlin/w0,modelfunCkLRZ(betaLCkLRZ(knnL(jL),:),wlin),'-',...
            'Color',cc(jL,:),'LineWidth',LineWidth)
    end
    set(gca,'FontSize',FontSize,'LineWidth',LineWidth,...
        'XLim',[0.1,0.7],...
        'XMinorTick','on','box','on')
    xlabel('$\omega$','Interpret','latex','FontSize',FontSize)
    ylabel('$C_L(k,\omega)$','Interpret','latex','FontSize',FontSize)
    
    %% Part 3
    axes('Position',[0.71 0.15 0.25 0.8])

    hold on
    plot(kbins,abs(betaLCkLRZ(:,1))/w0,'o-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','k',...
        'LineWidth',LineWidth)
    plot(kbins,abs(betaLCkLRZ(:,2))*pi/w0,'s-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','k',...
        'LineWidth',LineWidth)
    plot(kbins,abs(betaTCkLRZ(:,1))/w0,'o-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','w',...
        'LineWidth',LineWidth)
    plot(kbins,abs(betaTCkLRZ(:,2))*pi/w0,'s-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','w',...
        'LineWidth',LineWidth)
    patch([0,2,2,0],[wBP(1),wBP(1),wBP(2),wBP(2)],'k','FaceAlpha',.3)
    hold off
    box on
    set(gca,'XLim',[0,1.2],'YLim',[0,0.6],...
        'XMinorTick','on','FontSize',FontSize,'LineWidth',LineWidth)
    grid off
    xlabel('$k$','FontSize',FontSize,'Interpret','latex')
    ylabel('$\Omega(k), \pi\Gamma(k)$','FontSize',FontSize,'Interpret','latex')
    legend({'$\Omega_L$','$\pi\Gamma_L$','$\Omega_T$','$\pi\Gamma_T$'},...
        'Interpreter','latex','box','off','FontSize',FontSize,'Position',[0.9,0.25,0.01,0.01]);
    
    
    saveas(gcf,[path2.Results,filesep,'Fig3.fig'])
    saveas(gcf,[path2.Results,filesep,'Fig3.jpg'])
    saveas(gcf,[path2.Results,filesep,'Fig3.pdf'])
    saveas(gcf,[path2.Results,filesep,'Fig3.eps'])
end
toc