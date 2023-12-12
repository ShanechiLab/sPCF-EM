function [sigma,D,eigVecs] = getRandCov(dimXt,param1,param2,paramForm,diagFlag,eigVals)
    
    if ~exist('diagFlag','var'); diagFlag = false; end
    if ~exist('eigVals','var') || isempty(eigVals)
        if ~exist('paramForm','var'); paramForm = 'norm'; end
        switch paramForm
            case 'norm'
                center = param1;
                spread = param2;
                eigVals = abs(normrnd(center,spread,dimXt,1));
            case 'unif'
                if param1 > param2
                    error('wrong cov. matrix settings for uniformly generated')
                else
                    eigVals = getUnifRand(param1,param2,dimXt,1);
                end
        end
    end

    D = diag(eigVals);
    if diagFlag == true
        eigVecs = eye(length(eigVals));
    else
        eigVecs = orthCopy(randn(length(eigVals)));
    end
    sigma = zeros(dimXt,dimXt);
    for d = 1:dimXt
        v = eigVecs(:,d);
        sigma = sigma + (eigVals(d)*v)*v';
    end
    sigma = 0.5*(sigma+sigma');
end

