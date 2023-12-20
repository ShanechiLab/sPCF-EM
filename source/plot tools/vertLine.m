function h = vertLine(x,ax,LineType)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    if ~exist('ax','var') || isempty(ax)
        ax = gca;
    end
    if ~exist('LineType','var')
        LineType = 'k--';
    end
    
    hold on;
    curXlim = ax.XLim;
    curYlim = ax.YLim;
    curYRng = abs(diff(curYlim));
    h = plot([x,x],[curYlim(1)-3*curYRng, curYlim(2)+3*curYRng],LineType);
    ax.XLim = curXlim;
    ax.YLim = curYlim;
end

