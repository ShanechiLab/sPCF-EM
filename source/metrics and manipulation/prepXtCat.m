
function xCat = prepXtCat(anaS)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    if ~isfield(anaS,'stateXt')
        error('need xt to get xCat');
    end
    
    obsFields = {'obsNt'};
    found = false;
    ind = 1;
    nF = length(obsFields);
    while ind <= nF && ~found
        field = obsFields{ind};
        if isfield(anaS,field)
            found = true;
        else
            ind = ind + 1;
        end
    end
    if found
        obsField = field;
    else
        error('need observation fields');
    end
    
    if iscell(anaS.stateXt)
        dimXt = size(anaS.stateXt{1},1);
        nSegs = length(anaS.stateXt);
        xlen = 0;
        for s = 1:nSegs
            xSegLen = size(anaS.stateXt{s},2);
            obsSegLen = size(anaS.(obsField){s},2);
            if xSegLen == obsSegLen
                sInd = 1;
            else
                sInd = 2;
            end
            xlen = xlen + xSegLen - sInd + 1;
        end
        xCat = zeros(dimXt,xlen);
        ind = 1;
        for s = 1:nSegs
            xSegLen = size(anaS.stateXt{s},2);
            obsSegLen = size(anaS.(obsField){s},2);
            if xSegLen == obsSegLen
                sInd = 1;
            else
                sInd = 2;
            end
            indLen = xSegLen - sInd + 1;
            xCat(:,ind:ind+indLen-1) = ...
                            anaS.stateXt{s}(:,sInd:end);
            ind = ind + indLen;
        end
    else
        if size(anaS.stateXt,2)==size(anaS.(obsField),2)
            sInd = 1;
        else
            sInd = 2;
        end
        xCat = anaS.stateXt(:,sInd:end);
    end
end