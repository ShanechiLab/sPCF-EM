function out = centWithUpper(in,upper)
% CENTWITHUPPER Center input vector to 0 with upperbound on max
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    maxIn = max(in);
    meanIn = mean(in);
    
    if (0 - meanIn) > (upper - maxIn)
        out = (in - maxIn) + upper;
    else
        out = in - meanIn;
    end
end