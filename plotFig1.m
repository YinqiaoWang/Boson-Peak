dbstop if error
clearvars
close all

tic

%% Path
m = 2;

pathHome = 'E:\BosonPeak\BP4_1014';
folderList = {'180408a','180416a','180420a','180428a'};
nstepList = [10,7,6,1];

for jFolder = 1
    folderName = folderList{jFolder};
    
    path1 = [];
    folderListKt = [];
    load([pathHome,filesep,folderName,filesep,'path1.mat'],'path1','folderListKt');
    
    path2 = [];
    load([path1.data,filesep,folderListKt{m},filesep,'path2.mat'],'path2');
    
    nstep = nstepList(jFolder);
    for ksize = 1
        load([path2.dynamicMatrix,filesep,sprintf('%06d',nstep),'.',...
            sprintf('%02d',ksize),'.dynamicMatrix.mat'],'eigenfrequency')
        load([path2.dynamicStructureFactor,filesep,sprintf('%06d',nstep),'.',...
            sprintf('%02d',ksize),'.CurrentFunction2.mat'],'CLmean','CTmean','wbins','kbins','K')
        load([path2.modulus,filesep,sprintf('%06d',nstep),'.',sprintf('%02d',ksize),'.modulus.mat'],...
            'Gp_global','Gs_global','B_global','rho','vL','vT','diskNum','area')
        load([path2.modes,filesep,'ratioLow2.mat'],'ratioLow')
        load([path2.dynamicStructureFactor,filesep,sprintf('%06d',nstep),'.',...
            sprintf('%02d',ksize),'.SCkwFitting2.mat'],'wbins','betaTCwLRZ','betaLCwLRZ')
        
        
        w0 = 2220.8;
        
        dw = 100;   % \delta \omega
        nbins = 100; % number of bins
        wbinEdges = dw*(0:nbins);
        wbins0 = (wbinEdges(1:end-1)+wbinEdges(2:end))/2;
        
        eigenfrequency = eigenfrequency{1};
        
        wCounts = histcounts(eigenfrequency,wbinEdges);
        DOS0 = wCounts'/(dw*numel(eigenfrequency));
        DOSLow = cellfun(@sum,mat2cell(ratioLow,wCounts))/(dw*numel(eigenfrequency));
        DOSHigh = cellfun(@sum,mat2cell(1-ratioLow,wCounts))/(dw*numel(eigenfrequency));
        
%         kD = sqrt(4*pi*diskNum/area);
        kD = 3.52*10000/120;
        
        vT = wbins'./(abs(betaTCwLRZ(:,1))*10000/120);
        vL = wbins'./(abs(betaLCwLRZ(:,1))*10000/120);
        
        w = wbins(5:50);
        kT = abs(betaTCwLRZ(5:50,1))*10000/120;
        kL = abs(betaLCwLRZ(5:50,1))*10000/120;
        
        fitFunc = fittype('A*x/(w*sqrt(pi/2)).*exp(-2*(x-b).^2/w^2)+Q/pi*asin(pi*x/(Q*v))',...
            'independent',{'x'},'coefficients',{'A','w','b','Q','v'});
        kTw = fit(w',kT,fitFunc,'StartPoint', [5, 300, 900, 3.5*10000/120, 13]);
        kLw = fit(w',kL,fitFunc,'StartPoint', [5, 300, 900, 3.5*10000/120, 26]);
        
        wlin = (1:1000)';
        RDOS_pw = 1/kD^2*(kTw(wlin)./wlin.*differentiate(kTw,wlin)+kLw(wlin)./wlin.*differentiate(kLw,wlin));
        
    end
end

wbins = wbins0';
idxLow = wbins < 1500;

RDOS = DOS0./wbins;

wBP0 = 920;
fitRangew = find(wbins > 300 & wbins < 1500);
% wlin = linspace(wbins(fitRangew(1)),wbins(fitRangew(end)),1000);

modelfunLorentz = @(b,w) b(3)./(4*(w-b(1)).^2+b(2)^2);
betaLRZ0 = [wBP0,wBP0/2,max(RDOS)*wBP0^2];
betaLRZ = nlinfit(wbins(fitRangew),RDOS(fitRangew),modelfunLorentz,betaLRZ0);

wBP = (betaLRZ(1)+[-200,40])/w0;

plotScale = 2;

MarkerSize1 = plotScale*4;
MarkerSize2 = plotScale*3;
FontSize = plotScale*7;
LineWidth = plotScale*1;

cc = lines;

P = {[0.11,0.54,0.38,0.42];
    [0.11,0.12,0.38,0.42];
    [0.61,0.12,0.38,0.84]};

figure
set(gcf,'Units','centimeters','Position',[2,1,8,6]*plotScale)

ax1 = axes('Position',P{1});

% plot(wbins/w0,DOS0./wbins*w0^2,'o-','MarkerSize',MarkerSize1,...
%     'Color',cc(1,:),'MarkerFaceColor',cc(1,:),'LineWidth',LineWidth)
plot(ax1,wbins/w0,DOS0./wbins*w0^2,'-','LineWidth',LineWidth,'Color',cc(1,:))
patch([wBP(1),wBP(2),wBP(2),wBP(1)],[0,0,10,10],'k','FaceAlpha',.3,'EdgeColor','none')

set(ax1, 'TickDir', 'in','TickLength',[.02 .02], ...
    'XLim',[0,3],'YLim',[0,1.2],...
    'XMinorTick', 'on', 'YMinorTick', 'on',...
    'XTick',(0:1:3),'YTick',(0:0.4:1.2),...
    'XTickLabel',[],...
    'FontSize',FontSize,'LineWidth',LineWidth)

ax1.XAxis.MinorTickValues = (0.5:1:2.5);
ax1.YAxis.MinorTickValues = (0.2:0.4:1.0);

% xlabel(ax1,'$\omega$','FontSize',FontSize,'Interpret','latex')
ylabel(ax1,'$g(\omega)/\omega$','FontSize',FontSize,'Interpret','latex')

ax2 = axes('Position',P{2});

hold on
% plot(ax2,wbins/w0,DOS0*w0,'o','MarkerSize',MarkerSize2,'Color',cc(1,:),'MarkerFaceColor',cc(1,:))
plot(ax2,wbins/w0,DOS0*w0,'-','LineWidth',LineWidth,'Color',cc(1,:))
patch(ax2,[wBP(1),wBP(2),wBP(2),wBP(1)],[0,0,10,10],'k','FaceAlpha',.3,'EdgeColor','none')
hold off
box on
set(ax2, 'TickDir', 'in','TickLength',[.02 .02],...
    'XMinorTick', 'on', 'YMinorTick', 'on',...
    'XLim',[0,3],'YLim',[0,0.6],...
    'XTick',[0,1,2,3],'YTick',(0:0.2:0.4),...
    'FontSize',FontSize,'LineWidth',LineWidth)

ax2.XAxis.MinorTickValues = (0.5:1:2.5);
ax2.YAxis.MinorTickValues = (0.1:0.2:0.5);

xlabel(ax2,'$\omega$','FontSize',FontSize,'Interpret','latex')
ylabel(ax2,'$g(\omega)$','FontSize',FontSize,'Interpret','latex')

ax3 = axes('Position',P{3});
hold on
plot(ax3,wbins/w0,DOS0./wbins*w0^2,'o-','MarkerSize',MarkerSize1,...
    'Color',cc(1,:),'MarkerFaceColor',cc(1,:),'LineWidth',LineWidth)
plot(ax3,wbins(idxLow)/w0,DOSLow(idxLow)./wbins(idxLow)*w0^2,'^-','MarkerSize',MarkerSize1,...
    'Color',cc(2,:),'MarkerFaceColor',cc(2,:),'LineWidth',LineWidth)
plot(ax3,wbins(idxLow)/w0,DOSHigh(idxLow)./wbins(idxLow)*w0^2,'v-','MarkerSize',MarkerSize1,...
    'Color',cc(3,:),'MarkerFaceColor',cc(3,:),'LineWidth',LineWidth)
plot(ax3,wlin/w0,RDOS_pw*w0^2,'--','MarkerSize',MarkerSize1,...
    'Color','k','MarkerFaceColor','k','LineWidth',LineWidth)
patch(ax3,[wBP(1),wBP(2),wBP(2),wBP(1)],[0,0,10,10],'k','FaceAlpha',.3,'EdgeColor','none')
hold off

box on
set(ax3, 'TickDir', 'in','TickLength',[.02 .02],...
    'XMinorTick', 'on', 'YMinorTick', 'on',...
    'XLim',[0,0.7],'YLim',[0,1.1],...
    'XTick',[0,0.2,0.4,0.6],'YTick',(0:0.4:1.2),...
    'FontSize',FontSize,'LineWidth',LineWidth)

ax3.XAxis.MinorTickValues = (0.1:0.2:0.7);
ax3.YAxis.MinorTickValues = (0.2:0.4:1.0);

xlabel(ax3,'$\omega$','FontSize',FontSize,'Interpret','latex')
ylabel(ax3,'$g(\omega)/\omega$','FontSize',FontSize,'Interpret','latex')

legend({'$g(\omega)$','$g_a(\omega)$','$g_{na}(\omega)$','$g_{pw}(\omega)$'},...
    'Interpreter','latex','FontSize',FontSize,'box','off',...
    'Position',[0.34,0.75,0.1,0.15])

saveas(gcf,[path2.Results,filesep,'Fig1.fig'])
saveas(gcf,[path2.Results,filesep,'Fig1.jpg'])
saveas(gcf,[path2.Results,filesep,'Fig1.pdf'])
saveas(gcf,[path2.Results,filesep,'Fig1.eps'])



% figure
% set(gcf,'Units','centimeters','Position',[2,1,8,8]*plotScale)
%
% ax1 = axes('Position',[0.13,0.13,0.84,0.84]);
%
% set(ax1, 'TickDir', 'in','TickLength',[.02 .02], ...
%     'XLim',[0,3],'YLim',[0,1.2],...
%     'XTick',(0:0.5:3),'YTick',(0:0.2:1),...
%     'FontSize',FontSize,'LineWidth',LineWidth)
%
% box off
% ax2 = axes('Position', get(ax1, 'Position'), 'Color','none','LineWidth', LineWidth,...
%            'XAxisLocation','top', 'XTick', [],...
%            'YAxisLocation','right', 'YTick', []);
% linkaxes([ax1, ax2])
%
% ax2.NextPlot = 'add';
% plot(ax2,wbins/w0,DOS0./wbins*w0^2,'o-','MarkerSize',MarkerSize1,...
%     'Color',cc(1,:),'MarkerFaceColor',cc(1,:),'LineWidth',LineWidth)
% % plot(ax2,wlin/w0,modelfunLorentz(betaLRZ,wlin)*w0^2,'k--','LineWidth',LineWidth)
% plot(ax2,wbins(idxLow)/w0,DOSLow(idxLow)./wbins(idxLow)*w0^2,'^-','MarkerSize',MarkerSize1,...
%     'Color',cc(2,:),'MarkerFaceColor',cc(2,:),'LineWidth',LineWidth)
% plot(ax2,wbins(idxLow)/w0,DOSHigh(idxLow)./wbins(idxLow)*w0^2,'v-','MarkerSize',MarkerSize1,...
%     'Color',cc(3,:),'MarkerFaceColor',cc(3,:),'LineWidth',LineWidth)
% patch(ax2,[wBP(1),wBP(2),wBP(2),wBP(1)],[0,0,10,10],'k','FaceAlpha',.3,'EdgeColor','none')
%
% xlabel(ax1,'$\omega$','FontSize',FontSize,'Interpret','latex')
% ylabel(ax1,'$g(\omega)/\omega$','FontSize',FontSize,'Interpret','latex')
%
% legend({'$g(\omega)$','$g_a(\omega)$','$g_{na}(\omega)$'},...
%     'Interpreter','latex','FontSize',FontSize,...
%     'Position',[0.81,0.30,0.1,0.15])
%
% ax3 = axes('Position',[0.55,0.55,0.40,0.40]);
% set(ax3, 'TickDir', 'in','TickLength',[.02 .02],...
%     'XMinorTick', 'on', 'YMinorTick', 'on',...
%     'XLim',[0,3],'YLim',[0,0.6],...
%     'XTick',[0,1,2,3],'YTick',(0:0.2:0.6),...
%     'FontSize',FontSize,'LineWidth',LineWidth)
%
% box off
% ax4 = axes('Position', get(ax3, 'Position'), 'Color','none','LineWidth', LineWidth,...
%            'XAxisLocation','top', 'XTick', [],...
%            'YAxisLocation','right', 'YTick', []);
% linkaxes([ax3, ax4])
%
% ax4.NextPlot = 'add';
% plot(ax4,wbins/w0,DOS0*w0,'o','MarkerSize',MarkerSize2,'Color',cc(1,:),'MarkerFaceColor',cc(1,:))
% % plot(ax4,wbins(idxLow)/w0,DOSLow(idxLow)*w0,'^','MarkerSize',MarkerSize2,'Color',cc(2,:),'MarkerFaceColor',cc(2,:))
% % plot(ax4,wbins(idxLow)/w0,DOSHigh(idxLow)*w0,'v','MarkerSize',MarkerSize2,'Color',cc(3,:),'MarkerFaceColor',cc(3,:))
% patch(ax4,[wBP(1),wBP(2),wBP(2),wBP(1)],[0,0,10,10],'k','FaceAlpha',.3,'EdgeColor','none')
% % plot(ax3,wbins/w0,DOS2*w0./(wbins/w0),'-','Color','k','LineWidth',LineWidth)
% % plot(ax3,wbins/w0,DOST2*w0./(wbins/w0),'-','Color',cc(3,:),'LineWidth',LineWidth)
% % plot(ax3,wbins/w0,DOSL2*w0./(wbins/w0),'-','Color',cc(4,:),'LineWidth',LineWidth)
%
% % plot(ax3,[0,3,3],[0.6,0.6,0],'k-','LineWidth',LineWidth)
%
% % legend({'$g_0(\omega)$','$g(\omega)$','$g_T(\omega)$','$g_L(\omega)$'},...
% %     'Interpreter','latex','FontSize',FontSize,'box','off',...
% %     'Position',[0.82,0.3,0.1,0.1])
%
% set(ax3, 'XLim',[0,3],'YLim',[0,0.6],'LineWidth',LineWidth)
%
% ax3.XAxis.MinorTickValues = (0.5:1:2.5);
% ax3.YAxis.MinorTickValues = (0.1:0.2:0.5);
%
% xlabel(ax3,'$\omega$','FontSize',FontSize,'Interpret','latex')
% ylabel(ax3,'$g(\omega)$','FontSize',FontSize,'Interpret','latex')
%
%
% saveas(gcf,[path2.Results,filesep,'Fig1.fig'])
% saveas(gcf,[path2.Results,filesep,'Fig1.jpg'])
% % saveas(gcf,[path2.Results,filesep,'Fig1.pdf'])
% saveas(gcf,[path2.Results,filesep,'Fig1.eps'])
%
% figure
% set(gcf,'Units','centimeters','Position',[2,1,8,8]*plotScale)
% hold on
% plot(wbins/w0,DOS0*w0,'o','MarkerSize',MarkerSize2,'Color',cc(1,:),'MarkerFaceColor',cc(1,:))
% plot((wbins(idxLow)-gD*wbins(idxLow))/w0,DOSLow(idxLow)*w0,'^','MarkerSize',MarkerSize2,'Color',cc(2,:),'MarkerFaceColor',cc(2,:))
% plot(wbins(idxLow)/w0,DOSHigh(idxLow)*w0,'v','MarkerSize',MarkerSize2,'Color',cc(3,:),'MarkerFaceColor',cc(3,:))
% hold off
%
% set(gca,'XScale','log','YScale','log')

toc

