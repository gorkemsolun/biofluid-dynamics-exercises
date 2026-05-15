function plot_network(r,p,c)

        %% nice plot of hex network

        
        clf
        set(gcf,'color','white')
        set(gcf,'PaperPositionMode','auto')
        set(gcf,'InvertHardcopy','off')
        set(gca,'Position',[0.03 0.03 .95 .95])
        hold on
        plot([p(c(:,1),1),p(c(:,2),1)]',[p(c(:,1),2) p(c(:,2),2)]','k-','linewidth',16)
        plot([p(c(:,1),1),p(c(:,2),1)]',[p(c(:,1),2) p(c(:,2),2)]','w-','linewidth',10)
        scatter(r(:,1),r(:,2),100,'r.','sizedata',600);
        greenrbc = find(r(:,5)>0);
        scatter(r(greenrbc,1),r(greenrbc,2),100,'b.','sizedata',600);
        
        exitNode = find(p(:,1)+p(:,2) == max(p(:,1)+p(:,2)));
        scatter(p(exitNode,1),p(exitNode,2),100,'v','k');
        
        set(gca,'xticklabel','')
        set(gca,'color','white')
        set(gca,'yticklabel','')
        set(gca,'xcolor','white')
        set(gca,'ycolor','white')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'Ydir','reverse')
        axis([-0.03 1.03 -0.03 1.15]*max(p(:,1)))
        caxis([0 1])
        caxis([0 3e-4])
        axis equal
        hold off

        



drawnow
end