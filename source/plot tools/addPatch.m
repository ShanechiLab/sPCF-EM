function h = addPatch(ax,xs,ys,Color,alpha)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    if ~exist('ax','var') || isempty(ax)
        ax = gca;
    end
    if ~exist('xs','var') || isempty(xs)
        xs = ax.XLim;
    end
    if ~exist('ys','var') || isempty(ys)
        ys = ax.YLim;
    end
    
    xs = sort(xs);
    ys = sort(ys);
    
    Xdata = [xs(1),xs(1),xs(2),xs(2)];
    Ydata = [ys(1),ys(2),ys(2),ys(1)];
    
    if ~exist('alpha','var') || isempty(alpha)
        alpha = 1;
    end
    
    if ~exist('Color','var') || isempty(Color)
        h = patch('XData',Xdata,'Ydata',Ydata,'FaceAlpha',alpha,...
                                          'LineStyle','none');
    else
        h = patch('XData',Xdata,'Ydata',Ydata,'FaceColor',Color,...
                          'FaceAlpha',alpha,'LineStyle','none');
    end
    
end