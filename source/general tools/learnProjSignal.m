
function simTran = learnProjSignal(xTo,xFrom)
% Fit similarity transfrom from xFrom onto xTo
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    if size(xTo,2) ~= size(xFrom,2)
        error('must be same size');
    end
    
    toFrT = xTo*xFrom';
    frFrT = xFrom * xFrom';
    frFrTi = frFrT^-1;
    simTran = toFrT * frFrTi;
end