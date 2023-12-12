function stOut = convertSt(st,dimSt)
% CONVERTST convert markov chain to/from one-hot representations
%  Either converts st from a 1xT array of values of [1:M] to an equivalent
%  MxT matrix of one-hot vectors, or vice-versa.
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    if size(st,1) > size(st,2)
        st = st';
        flipAfter = true;
    else
        flipAfter = false;
    end
    
    if size(st,1) == 1
        if ~exist('dimSt','var')
            dimSt = length(unique(st));
        end
        
        if dimSt > 1
            stOut = zeros(dimSt,size(st,2));
            for k = 1:dimSt
                stOut(k,st == k) = 1;
            end
        else
            stOut = st;
        end
    else
        stOut = zeros(1,size(st,2));
        for k = 1:size(st,1)
            stOut(logical(st(k,:))) = k;
        end
    end
    if flipAfter == true
        stOut = stOut';
    end
end

