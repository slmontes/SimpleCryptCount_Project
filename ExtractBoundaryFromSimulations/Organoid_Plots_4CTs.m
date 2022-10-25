%Plot the epithelial layer of an organoid over time

function [epithelialCells, nonGhostsx, nonGhostsy,panethCells,stemCells,TACells,ECcells] = Organoid_Plots_4CTs(ghostx,ghosty,tempx,tempy,panethx,panethy,stemx,stemy,TAx,TAy,ECx,ECy,timePoint)
%%
% NEED: ghostx,ghosty,tempx,tempy,timePoint
%
% CanGIVE: epithelialCells, nonGhostsx, nonGhostsy. panethCells, stemCells,
% TACells, EC cells
%%
epithelialCells = [tempx(timePoint,:)',tempy(timePoint,:)'];
epithelialCells(~any(epithelialCells,2),:)=[];

nonGhostsx = [ghostx(timePoint,:)';epithelialCells(:,1)];%all matrigel and organoid cells
nonGhostsy = [ghosty(timePoint,:)';epithelialCells(:,2)];

%The following figures show a plot of all non-ghost cells found in the
%simulation as empty black circles and the epithelial cells as magenta stars.
figure('visible','off');hold on;
plot(nonGhostsx,nonGhostsy,'ko');
% plot(epithelialCells(:,1),epithelialCells(:,2),'m*');
%%
%To avoid errors in case there is any of these cell types missing we have
%to define them:

panethCells = [panethx(timePoint,:)',panethy(timePoint,:)'];
stemCells = [stemx(timePoint,:)',stemy(timePoint,:)'];
TACells = [TAx(timePoint,:)',TAy(timePoint,:)'];
ECcells = [ECx(timePoint,:)',ECy(timePoint,:)'];

panethCells(~any(panethCells,2),:)=[]; %remove all (0,0) points from organoid
stemCells(~any(stemCells,2),:)=[];
TACells(~any(TACells,2),:)=[];
ECcells(~any(ECcells,2),:)=[];

%The following plot shows not only the epithelial layer but each cell type
%that forms it
figure('visible','off'); hold on;
plot(stemCells(:,1),stemCells(:,2),'b*');
plot(panethCells(:,1),panethCells(:,2),'g*');
plot(TACells(:,1),TACells(:,2),'y*');
plot(ECcells(:,1),ECcells(:,2),'*','MarkerFaceColor',[.5 0 .5]);
axis([0 30 0 30]);
title('2D Intestinal Organoid In-silico');
legend('Stem cells','Paneth cells','TA cells','EC cells')
set(gcf,'Position',  [500, 500, 1000, 800]);
set(gca, 'FontSize',13, 'XTick',[], 'YTick',[]);

end
