function [Gmag,Gdir,Gx,Gy]=plotGmag(img,tp3new,tempid,plotid,option,s5v2)

switch option
    case 'no plot'
        [Gx, Gy] = imgradientxy(img);
        [Gmag, Gdir] = imgradient(Gx, Gy);        
        
    case 'plot'
        
        [Gx, Gy] = imgradientxy(img);
        [Gmag, Gdir] = imgradient(Gx, Gy);
        
        x=diff(s5v2);
        x = x(tempid);
        x = scaledata(x,min(sum(Gmag)), max(sum(Gmag)));
        
        subplot(2,length(tp3new),plotid);imshowpair(img,Gmag,'montage'); axis off;
        
        subplot(2,length(tp3new),plotid+length(tp3new));
        plot(tempid,sum(Gmag),'LineWidth',2);hold on;plot(tempid,sum(img),'m');hold on;plot(tempid,x,'k');
        xl = get(gca,'XLim');
        if length(tempid)<=6
            set(gca,'XTick',linspace(xl(1),xl(2),length(tempid)),'XTickLabel',tempid)
        else
            set(gca,'XTickLabel',tempid)
        end
        
        clear xl
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        
        text(0.5, 1,'\bf Original image (left); Gradient magnitude (right)','HorizontalAlignment' ,'center','VerticalAlignment', 'top')
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
end