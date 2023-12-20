function newModel = trimModel(mdlS,finInd)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    nIter = length(mdlS.thetaCell)-1;
    if ~exist('finInd','var')
        finInd = getFinalCellInd(mdlS.thetaCell) - 1;
    end
    
    if finInd < nIter
        fields = setdiff(fieldnames(mdlS),{'thetaCell','static'});
        newModel = struct;
        newModel.thetaCell = mdlS.thetaCell(1:finInd+1);
        newModel.static = mdlS.static;
        for fInd = 1:length(fields)
            field = fields{fInd};
            newModel.(field) = mdlS.(field)(:,1:finInd);
        end
    else
        newModel = mdlS;
    end
end