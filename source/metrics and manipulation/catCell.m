
function xCat = catCell(xC)

    if iscell(xC)
        nC = length(xC);
        totLen = 0;
        for c = 1:nC
            totLen = totLen + size(xC{c},2);
        end
        ord = size(xC{1},1);
        xCat = NaN(ord,totLen);
        ind = 0;
        for c = 1:nC
            len = size(xC{c},2);
            xCat(:,ind+1:ind+len) = xC{c};
            ind = ind + len;
        end
    else
        xCat = xC;
    end
end