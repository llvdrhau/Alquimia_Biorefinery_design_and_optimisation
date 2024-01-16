function []=f_my_contorno(matrix,pH,HRT)
% close all
% figure


    
A=matrix;
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
% fig1 = figure;
Y = pH;
X = HRT;
[c,h1]=contourf(X,Y,A);
% caxis([0 7])
% clabel(c,h1,'FontSize',12);
% hLabel=clabel(c, h1);
% set(hLabel,'FontName','Garamond','FontSize',20,'Rotation',0,'Color','r')
colorbar;
%caxis([0,50])
hold on
set(gca,'Color','k')
ax.FontSize = 2;
% ax.Color = 'r';
xName = {'HRT (days)'};
% yName = 'Solids exchange ratio';
yName = 'pH';
xlabel(xName,'Color','k')
ylabel(yName,'Color','k')
grid off
%title('Butyrate concentration (C-mol/C-mol feeding)','FontSize',16)
% set(gca,'xtick',[0 20 40 60 80 100])
ylim([5 8])
set(gca,'ytick',[[5:1:8]])
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))
xlim([2 6])
set(gca,'xtick',[2:6])
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.0f'))


ax1 = gca;
ax1.FontSize = 16;
ax1.Color = 'k';
% set(gcf,'position',[10,10,1400,900])
ax1_pos = ax1.Position;
ax1.TickDir = 'out';

% textoA = 'Biomass washout';
% h=text(4,0.19,textoA,'FontSize',32,'color','white');
% set(h,'Rotation',73)
% % % % 
% textoB = 'CH_4';
% text(20,0.12,textoB,'FontSize',22,'color','white')
% % % 
% textoC = 'formation';
% text(17.7,0.11,textoC,'FontSize',22,'color','white')
% 
% ax2 = axes('Position',ax1_pos,'XAxisLocation','top','Color','none');
% ax2 = axes('XAxisLocation','top','Color','none');
% set(gca,'ytick',[])
% 
% xlabel('Glucose concentration (g/L)','FontSize',10,'Color','k')
% ax2.TickDir = 'out';
% ax2.FontSize = 20;
% set(gca,'ytick',[])
% xlim([0 100])
% set(gca,'xtick',[X])
% set ( gca, 'xdir', 'reverse' )
%      grid off
% %  set(gcf, 'PaperPosition', [0 0 20 20])    % can be bigger than screen 
% % set(gcf, 'PaperSize', [20 20])
% set(gcf,'InvertHardCopy','off')
% grid off

% print(gcf, 'nBut2.jpg', '-djpeg', '-r300' );