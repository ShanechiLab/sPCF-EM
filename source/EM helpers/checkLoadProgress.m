function [complete,loadFinIter,newIter] = checkLoadProgress(loadMdl,nIter)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    ldStc = loadMdl.static;
    
    complete = false;
    loadSetIter = length(loadMdl.thetaCell)-1;
    loadFinIter = getFinalCellInd(loadMdl.thetaCell) - 1;
    newIter = [];

    if isfield(ldStc,'bailedOut')
        complete = true;
    else
        if nIter <= loadSetIter
            complete = ~isempty(loadMdl.thetaCell{nIter+1});
        else
            newIter = nIter;
        end
    end
end