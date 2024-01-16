function []=f_my_contorno2(matrix,substrates_ratio,pH)
    
A=matrix;
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

Y = pH;
X =substrates_ratio*100;
[c,h1]=contourf(X,Y,A);

colorbar;

hold on
set(gca,'Color','k')
ax.FontSize = 16;
xName = {'TWW (% COD)'};
yName = 'pH';
xlabel(xName,'Color','k')
ylabel(yName,'Color','k')
grid off

ylim([5 8])
set(gca,'ytick',[[5:1:8]])
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))
xlim([0 100])
set(gca,'xtick',[0:20:100])
tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.0f'))

ax1 = gca;
ax1.FontSize = 16;
ax1.Color = 'k';
ax1_pos = ax1.Position;
ax1.TickDir = 'out';

ax2 = axes('Position',ax1_pos,'XAxisLocation','top','Color','none');
ax2.YTick=[];
ax2.FontSize = 16;
linkprop([ax1, ax2],{'Units','Position','ActivePositionProperty'});
xlim([0 100]);
set(gca,'xtick',[0:20:100]);
ax2.XDir='reverse';
xlabel('Regrind pasta (% COD)','Color','k')

set(gcf,'position',[5,5,800,600])


end
