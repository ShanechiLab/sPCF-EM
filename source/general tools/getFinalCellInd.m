function finalInd = getFinalCellInd(myCell)
% Get index of final non-empty element in cell array
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    notValid = ~iscell(myCell) || length(size(myCell)) > 2 ...
                               || min(size(myCell)) > 1;
    if ~notValid
        totLen = length(myCell);
        firstEmptyInd = find(cellfun('isempty',myCell),1);
        if isempty(firstEmptyInd)
            finalInd = totLen;
        else
            finalInd = firstEmptyInd-1;
        end
    else
        finalInd = NaN;
    end
end