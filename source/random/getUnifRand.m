function out = getUnifRand(low,high,varargin)

    if isempty(high) && length(low) == 2
        high = max(low);
        low = min(low);
    end

    range = abs(low-high);
    low = min(low,high);
    if isempty(varargin)
        out = range*rand(1) + low;
    else
        out = range*rand([varargin{:}]) + low;
    end
end