function out = isfinitereal(A)
    out = all(reshape(isfinite(A) & isreal(A),1,[]));
    
end