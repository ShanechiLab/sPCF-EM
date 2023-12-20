function PNew = rescaleDist(PStEst,factor)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    [dimSt,tlen] = size(PStEst);
    
    if ~exist('factor','var')
        factor = 1;
    end
    
    if factor ~= 1
        tmp = PStEst.^factor;
        PNew = (1./sum(tmp,1)).*tmp;
    else
        pNew = PStEst;
    end

end