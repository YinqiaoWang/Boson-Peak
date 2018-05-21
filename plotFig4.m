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
    
    load([path2.modulus,filesep,sprintf('%06d',nstep),'.',sprintf('%02d',ksize),'.modulus.mat'],...
        'Gp_global','Gs_global','B_global','rho','vL','vT');
    
    load([path2.dynamicStructureFactor,filesep,sprintf('%06d',nstep),'.',...
        sprintf('%02d',ksize),'.SCkwFitting2.mat'],...
        'wbins','kbins','CL','CT',...
        'modelfunCkLRZ','betaTCkLRZ','betaLCkLRZ',...
        'betaTCwLRZ','betaLCwLRZ')
    
    
    w0 = 2220.8;
    wBP = [720,1120]/w0;
    
    vL0 = vL*10000/(51+70);    % m/s -> D/s
    vT0 = vT*10000/(51+70);    % m/s -> D/s
    
    dw = 100;   % \delta \omega
    nbins = 70; % number of bins
    wbinEdges = dw*(0:nbins);
    wbins0 = (wbinEdges(1:end-1)+wbinEdges(2:end))/2;
    
    plotScale = 2;
    MarkerSize = plotScale*4;
    FontSize = plotScale*7;
    LineWidth = plotScale*1;
    
    cc = lines;
    
    P = {[0.1,0.6,0.38,0.38];
         [0.595,0.6,0.38,0.38];
         [0.1,0.1,0.38,0.38];
         [0.595,0.1,0.38,0.38]};
    figure
    set(gcf,'Units','centimeters','Position',[2,2,8,8]*plotScale)
    
    %% Ioffe-Regel
    axes('Position',P{1})

    hold on
    plot(kbins,abs(betaLCkLRZ(:,1))/w0,'o-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','k',...
        'LineWidth',LineWidth)
    plot(kbins,abs(betaLCkLRZ(:,2))*pi/w0,'s-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','k',...
        'LineWidth',LineWidth)
    plot(kbins,abs(betaTCkLRZ(:,1))/w0,'o-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','w',...
        'LineWidth',LineWidth)
    plot(kbins,abs(betaTCkLRZ(:,2))*pi/w0,'s-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','w',...
        'LineWidth',LineWidth)
    patch([0,2,2,0],[wBP(1),wBP(1),wBP(2),wBP(2)],'k','FaceAlpha',.3,'EdgeColor','none')
    hold off
    box on
    set(gca,'XLim',[0,1.2],'YLim',[0,0.6],...
        'XMinorTick','on','FontSize',FontSize,'LineWidth',LineWidth)
    grid off
    xlabel('$k$','FontSize',FontSize,'Interpret','latex')
    ylabel('$\Omega(k), \pi\Gamma(k)$','FontSize',FontSize,'Interpret','latex')
    legend({'$\Omega_L$','$\pi\Gamma_L$','$\Omega_T$','$\pi\Gamma_T$'},...
        'Interpreter','latex','box','off','FontSize',FontSize,'Position',[0.406,0.67,0.01,0.01]);
    
    
    %% Velocity dip
    vT = wbins'./abs(betaTCwLRZ(:,1));
    vL = wbins'./abs(betaLCwLRZ(:,1));
    
    wbin_counts = histcounts(wbins,wbinEdges);
    vT = cellfun(@sum,mat2cell(vT,wbin_counts))*20/dw;
    vL = cellfun(@sum,mat2cell(vL,wbin_counts))*20/dw;
    
    ax1 = axes('Position',P{3});
    
    hold on
    plot(wbins0(3:end)/w0,vL(3:end)/vL0,'o-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','k',...
        'LineWidth',LineWidth)
    plot(wbins0(3:end)/w0,vT(3:end)/vT0,'s-','Color','k','MarkerSize',MarkerSize,'MarkerFaceColor','w',...
        'LineWidth',LineWidth)
    patch([wBP(1),wBP(1),wBP(2),wBP(2)],[0,10,10,0],'k','FaceAlpha',.3,'EdgeColor','none')
    hold off
    
    set(gca,'FontSize',FontSize,'box','on','LineWidth',LineWidth,...
        'XLim',[0,0.8],'YLim',[0.6,1],...
        'XTick',(0:0.2:0.8),'YTick',[0.6,0.8,1],...
        'XMinorTick','on','YMinorTick','on')
    ax1.XAxis.MinorTickValues = (0.1:0.2:0.7);
    ax1.YAxis.MinorTickValues = [0.7,0.9];
    
    xlabel('$\omega$','Interpret','latex')
    ylabel('$v(\omega)/v_0$','FontSize',FontSize,'Interpret','latex')
    legend({'$v_L$','$v_T$'},'FontSize',FontSize,'Interpreter','latex',...
        'box','off','Position',[0.41,0.42,0.01,0.01])
    

    %% Participation Ratio
    p = [];
    w = [];
    load([path2.modes,filesep,sprintf('%06d',nstep),'.',sprintf('%02d',ksize),...
        '.participationRatio.mat'],'p','w')
    
    pmean = zeros(nbins,1);
    
    for wn = 1:nbins
        idx_wn = w >= dw*(wn-1) & w < dw*wn;
        pmean(wn) = mean(p(idx_wn));
    end
    
    ax2 = axes('Position',P{2});
    
    hold on
    plot(w/w0,p,'.','Color',cc(1,:))
    plot(wbins0/w0,pmean,'-','Color',cc(2,:),'LineWidth',LineWidth)
    patch([wBP(1),wBP(1),wBP(2),wBP(2)],[1e-4,10,10,1e-4],'k','FaceAlpha',.3,'EdgeColor','none')
    hold off
    set(gca, 'Box', 'on','LineWidth',LineWidth,  ...
        'XMinorTick', 'on', ...
        'XTick',(0:0.2:0.8),...
        'XLim',[0,0.8],'YLim',[0,0.7],...
        'FontSize',FontSize)
    ax2.XAxis.MinorTickValues = (0.1:0.2:0.7);
    
    xlabel('$\omega$','Interpret','latex')
    ylabel('$p(\omega)$','Interpret','latex')
    
    
    
    %% Diffusivity
    w = [];
    diffusivity = [];
    load([path2.dynamicMatrix,filesep,sprintf('%06d',nstep),'.',...
        sprintf('%02d',ksize),'.diffucivity.mat'],'diffusivity','w')
    
    
    diffusivityMean = zeros(nbins,1);
    for jbins = 1:nbins
        idx_wbins = w > dw*(jbins-1) & w <= dw*jbins;
        diffusivityMean(jbins) = mean(diffusivity(idx_wbins));
    end
    
    ax3 = axes('Position',P{4});
    
    hold on
    plot(w/w0,diffusivity,'.','Color',cc(1,:))
    plot(wbins0/w0,diffusivityMean,'-','Color',cc(2,:),'LineWidth',LineWidth)
    patch([wBP(1),wBP(1),wBP(2),wBP(2)],[1e-4,10,10,1e-4],'k','FaceAlpha',.3,'EdgeColor','none')
    hold off
    set(gca,'box','on','LineWidth',LineWidth,...
        'XLim',[0,0.8],'YLim',[0,1],...
        'XTick',(0:0.2:0.8),'YTick',(0:0.4:0.8),...
        'XMinorTick','on','YMinorTick','on',...
        'FontSize',FontSize)
    ax3.XAxis.MinorTickValues = (0.1:0.2:0.7);
    ax3.YAxis.MinorTickValues = [0.2,0.6];
    xlabel('$\omega$','Interpret','latex')
    ylabel('$d(\omega)$','Interpret','latex')

    
    saveas(gcf,[path2.Results,filesep,'Fig4.fig'])
    saveas(gcf,[path2.Results,filesep,'Fig4.jpg'])
    saveas(gcf,[path2.Results,filesep,'Fig4.pdf'])
    saveas(gcf,[path2.Results,filesep,'Fig4.eps'])
end
toc

