function finalInd = getFinalCellInd(myCell)
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