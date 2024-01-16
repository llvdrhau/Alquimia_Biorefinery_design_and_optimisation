function f_plot_my_results(tOut,y,compounds,colour,legendText)

% This function plots the VFA production results
% *********************************************************************** %
% INPUTS:
% time ~ tOut
% states ~ y
% compounds
% colour
% legendText
% *********************************************************************** %
% OUTPUT:
% figures with the plotted results
% *********************************************************************** %

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela.
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code. 


if nargin<4
    colour='b'; %By default blue
end

%Set up default plot caracteristics
set(0,'DefaultLineLineWidth',1,'DefaultAxesFontSize',11,...
    'DefaultAxesXGrid','off','DefaultAxesYGrid','off',...
    'DefaultLineMarkerSize',4.5,...
    'DefaultAxesFontName','Century Schoolbook',...
    'DefaultAxesXcolor', 0.25*[1,1,1],'DefaultAxesYcolor', 0.25*[1,1,1]);

states=compounds.Abb;
indexstates=compounds.index;

h =  findobj('type','figure');
nPreviousFig = length(h);
nPlots=length(indexstates);  %As many plots as states
nFigures=ceil(nPlots/4);        %Plots grouped in 4 subplots figures
iCompound=1;
for i=nPreviousFig+1:nFigures+nPreviousFig
    figure(i)
    for j=1:4
        if iCompound>nPlots     
            break
        end
        subplot(2,2,j)
        yVector=y(:,iCompound); %Select state
        plot(tOut,yVector,colour)
        hold on
        xlabel('Time (days)')
        yLabelText=horzcat(states{iCompound},'  (g/L)');
        ylabel(yLabelText)
        hold on
        iCompound=iCompound+1;
    end
    if nargin>4
        legend(legendText)
    end
    
end