function plot3d_optimization_results(allHPLCData,xlab,ylab,zlab,titl)
   n = length(allHPLCData(:,1));
   figure('Position', [100, -100, 1049, 895])
   %allHPLCData(:,7) = 1 - allHPLCData(:,7);
   maxV = 0.0224;%min(allHPLCData(:,7));
   %index = find(allHPLCData(:,7)==1);
   %allHPLCData(index,7) = 0;
   k = 100;
   %c = colormap(flipud(jet(k)));
   c = colormap(jet(k));
   y = linspace(.0224, .4);
   for l = 1:n
      ind = find(y>allHPLCData(l,7));
      ind = 100-ind(1);
      stem3(allHPLCData(l,3)/allHPLCData(l,4),allHPLCData(l,5), allHPLCData(l,2),'filled','Color',c(ind,:),'MarkerSize',15)
      hold on 
   end
   xlabel(xlab,'FontSize', 20);
   ylabel(ylab,'FontSize', 20);
   zlabel(zlab,'FontSize', 20);
   title(titl,'FontSize', 30);
   set(gca,'fontsize',14);
   colorbar('Ticks',[0,.5,1],...
         'TickLabels',{'Low','Mid','High'})
   set(gcf,'color','w');
   hold off
   
%      x=rand(20,1); %data to be plotted
% ran=range(x); %finding range of data
% min_val=min(x);%finding maximum value of data
% max_val=max(x)%finding minimum value of data
% y=floor(((x-min_val)/ran)*63)+1; 
% col=zeros(20,3)
% p=colormap(hsv(64))
% for i=1:20
%   a=y(i);
%   col(i,:)=p(a,:);
%   stem3(i,i,x(i),'Color',col(i,:))
%   hold on
% end
end