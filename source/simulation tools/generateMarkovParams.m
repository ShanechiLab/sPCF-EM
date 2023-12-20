function [sTran,sInit] = generateMarkovParams(pStay,dimSt)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    if dimSt <= 1
        if dimSt < 1
            error('dimSt must be greater than 0');
        else
            sInit = 1;
            sTran = 1;
        end
    else
        pSwitch = (1-pStay)/(dimSt-1);
        tRow = pSwitch.*ones(1,dimSt);
        tRow(1) = pStay;
        sTran = toeplitz(tRow);
        sInit = (1/dimSt).*ones(dimSt,1);
    end
    
    sTran = sTran + diag(1-sum(sTran,1));
    sInit(1) = sInit(1) + (1 - sum(sInit));
end