function [Norm_Aseg1, CryptSection, Curv_sum] = CalculateSectionArea(ax,ay,bx,by,x,y,k,N)
dbstop if error
% a and b are points obtained from the midrow_section vector
Line1x = [ax, bx]';
Line1y = [ay, by]';

IDXx = knnsearch(x,Line1x);
IDXy = knnsearch(y,Line1y);

% if IDXx~=IDXy
%     IDXy = IDXx;
% end

%Select entire crypt segment 
seg1x = x(IDXx(1):IDXx(2));
seg1y = y(IDXy(1):IDXy(2));

%Select the calculated curvature for the specific segment points 
Curv_seg = k(IDXx(1):IDXx(2));
Normals_kx = N((IDXx(1):IDXx(2)),1);
Normals_ky = N((IDXy(1):IDXy(2)),2);

if IDXx(1)>IDXx(2)
    seg1x = [x(IDXx(1):end);x(1:IDXx(2))];
    seg1y = [y(IDXy(1):end);y(1:IDXy(2))];
    
    Curv_seg = [k(IDXx(1):end);k(1:IDXx(2))];
    Normals_kx = [N((IDXx(1):end),1);N((1:IDXx(2)),2)]; 
    Normals_ky = [N((IDXy(1):end),1);N((1:IDXy(2)),2)];
end

CryptSection = [seg1x,seg1y];

Curv_sum = abs(sum(Curv_seg));

ATotal = polyarea(x,y);
Aseg1 = polyarea(seg1x,seg1y);

%Normalize area according the total organoid area
Norm_Aseg1=Aseg1/ATotal;

%% UNComment this senction to record a video!
% % figure
% hold on
% plot(x, y,'r-','LineWidth',0.5);
% % plot(seg1x, seg1y,'b-','LineWidth',0.5);
% plot(CryptSection(:,1), CryptSection(:,2),'b-','LineWidth',0.5);
% plot(Line1x, Line1y,'g')
% % plot([seg1x, seg1x+Curv_seg.*Normals_kx]', [seg1y, seg1y+Curv_seg.*Normals_ky]','m','LineWidth',2.5);
% % hold off
% % axis([-150 150 -150 150]);
% axis([-200 200 -150 150]);
% % axis([-300 300 -300 300]);
% title('2D Intestinal Organoid In-vitro Day 7 -Sectioning');
% % legend('Stem cells','Paneth cells','TA cells')
% set(gcf,'Position',  [500, 500, 1000, 800]);
% set(gca, 'FontSize',13, 'XTick',[], 'YTick',[]);

end