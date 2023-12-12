
function [PPavg,PPsep,AUCavg,AUCsep] = getPP(PNtAll,ntAll,field)

    if ~iscell(PNtAll)
        PNtAll = {PNtAll};
        ntAll = {ntAll};
    end
    nSegs = length(PNtAll);
    if isstruct(PNtAll{1}) && exist('field','var')
        pIn = PNtAll;
        for s = 1:nSegs
            PNtAll{s} = pIn{s}.(field);
        end
    end
    
    ntTot = catCell(ntAll);
    PNtTot = catCell(PNtAll);
    dimNt = size(ntTot,1);
    
    PPsep = NaN(dimNt,1);
    for c = 1:dimNt
        if sum(ntTot(c,:)) > 0
            auc = scoreAUC(logical(ntTot(c,:))',PNtTot(c,:)');
            PPsep(c) = 2*auc - 1;
        end
    end

    PPavg = nanmean(PPsep);
    AUCsep = (PPsep + 1)./2;
    AUCavg = nanmean(AUCsep);
end