function hitRate = calcStAccuracy(PStEst,stTrue)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    if size(stTrue,1) > size(stTrue,2)
        stTrue = stTrue';
    end

    if size(stTrue,1) > 1
        dimStTru = size(stTrue,1);
        stTrue = convertSt(stTrue);
    else
        dimStTru = max(stTrue);
    end
    % stTrue is array of true value of regime
    N = length(PStEst);
    hitRate = 0;
    if length(PStEst) == length(stTrue)
        if size(PStEst,1) > size(PStEst,2)
            PStEst = PStEst';
        end

        if size(PStEst,1) > 1
            
            dimStEst = size(PStEst,1);
            if dimStEst == dimStTru
                hitRate = zeros(factorial(dimStTru),1);
                allStPerms = perms(1:dimStTru);
                for i = 1:factorial(dimStTru)
                    stPerm = allStPerms(i,:);
                    [~,stEst] = max(PStEst(stPerm,:),[],1);
                    hitRate(i) = sum(stEst == stTrue) ./ N;
                end
            end
        else
            stEst = PStEst;
            hitRate = sum(stEst == stTrue) ./ N;
        end
    end
end