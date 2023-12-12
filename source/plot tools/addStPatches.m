function addStPatches(t,st,yvals,colors,alpha)

    if ~exist('alpha','var')
        alpha = 0.5;
    end

    ax = gca;
    tlen = length(t);
    edges = [1,find(abs(diff(st))),tlen];
    
    tEdges = t(edges);
    nPatches = length(edges) - 1;
    
    if ~exist('yvals','var') || isempty(yvals)
        yvals = ax.YLim;
    end
    
    for k = 1:nPatches
       tStart = tEdges(k);
       tEnd = tEdges(k+1);
       sVal = st(edges(k)+1);
       if sVal > 0
           col = colors(sVal,:);
           addPatch(ax,[tStart,tEnd],yvals,col,alpha);
       end
       hold on;
    end
    


end