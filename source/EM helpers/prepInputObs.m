function [obsNtAll,ntOut,dimNt] = prepInputObs(trnS,prmS,noCell)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    if ~exist('prmS','var')
        prmS = struct;
    end
    if ~exist('noCell','var')
        noCell = false;
    end
    
    varList = {};
    defaultVals = {};    

    prmS = insertStructDefaults(prmS,varList,defaultVals);
    
    ntIsCell = false;
    if isfield(trnS,'obsNt') && ~isempty(trnS.obsNt)
        if iscell(trnS.obsNt)
            finiteFound = false;
            ind = 1;
            while ind <= length(trnS.obsNt) && ~finiteFound
                if any(isfinite(trnS.obsNt{ind}(:)))
                    finiteFound = true;
                end
                ind = ind + 1;
            end
        else
            finiteFound = any(isfinite(trnS.obsNt(:)));
        end
        if finiteFound == true
            ntFound = true;
            if iscell(trnS.obsNt)
                ntIsCell = true;
            end
        else
            ntFound = false;
        end
    else
        ntFound = false;
    end

    if ~ntFound
        error('training struct must contain obsNt');
    end

    obsCell = ntIsCell;
    if obsCell == false
        ntOrig = {trnS.obsNt};

        dimNt = size(ntOrig{1},1);
    else
        ntOrig = trnS.obsNt;
        dimNt = size(ntOrig{1},1);
    end
    
    nSegs = length(ntOrig);
    obsNtAll = cell(1,nSegs);
    ntOut = ntFound;
    for sInd = 1:nSegs
        obsNtAll{sInd} = ntOrig{sInd};
    end
    if nSegs == 1 && noCell
        obsNtAll = obsNtAll{1};
    end
end
