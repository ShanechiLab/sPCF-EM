function newS = expandModel(mdlS,nIterNew)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    nIter = length(mdlS.thetaCell)-1;

    if nIter < nIterNew
        fields = setdiff(fieldnames(mdlS),{'thetaCell','static'});
        newS = struct;
        newS.thetaCell = cell(1,nIterNew+1);
        newS.thetaCell(1:nIter+1) = mdlS.thetaCell;
        newS.static = mdlS.static;
        for fInd = 1:length(fields)
            field = fields{fInd};
            fieldSize = size(mdlS.(field),1);
            newS.(field) = zeros(fieldSize,nIterNew);
            newS.(field)(:,1:nIter) = mdlS.(field);
        end
    else
        newS = mdlS;
    end
end
