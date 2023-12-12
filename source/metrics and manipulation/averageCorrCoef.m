function [cc,ccs] = averageCorrCoef(xt,zt)

    [dims,samps] = size(xt);
    
    ccs = zeros(dims,1);
    for k = 1:dims
        x = xt(k,:)';
        y = zt(k,:)';
        R = corrcoef(x,y);
        ccs(k) = R(1,2);
    end

    if any(isnan(ccs))
        ccs(isnan(ccs)) = 0;
    end
    cc = mean(ccs);


end