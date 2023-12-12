
function simTran = learnProjSignal(xTo,xFrom)

    if size(xTo,2) ~= size(xFrom,2)
        error('must be same size');
    end
    
    toFrT = xTo*xFrom';
    frFrT = xFrom * xFrom';
    frFrTi = frFrT^-1;
    simTran = toFrT * frFrTi;
end