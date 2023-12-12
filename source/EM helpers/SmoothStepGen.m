function [xEst,KEst,J] = SmoothStepGen(xDec,KDec,xDif,KDif,Jright,A,Q,...
                                       xSmt,KSmt,xPrd,KPrd,KPrdi)
% SMOOTHSTEPGEN Rauch-Tung-Striebel (RTS) Kalman smoother step
%  One backwards step of the RTS smoother. Different inputs possible based
%  on what has potentially already been calculated. General equation takes
%  the following form:
%       xEst = xDec + J*(xSmt-xPrd)
%       KEst = KDec + J*(KSmt-KPrd)*J'
%          J = KDec*A'*(KPrd^-1);
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    needDiff = ~exist('xDif','var') || isempty(xDif);
    noJright = ~exist('Jright','var') || isempty(Jright);
    noPrd = ~exist('xPrd','var') || isempty(xPrd);
    noPrdI = ~exist('KPrdi','var') || isempty(KPrdi);
    
    if noPrd
        if needDiff || (noJright && noPrdI)
            xPrd = A*xDec;
            KPrd = A*KDec*A' + Q;
        end
    end

    if needDiff
        xDif = xSmt - xPrd;
        KDif = KSmt - KPrd;
    end

    if noJright
        if noPrdI
            J = KDec*A'/KPrd;
        else
            J = KDec*A'*KPrdi;
        end
    else
        J = KDec*Jright;
    end

    xEst = xDec + J*xDif;
    Ktmp = KDec + J*KDif*J';
    KEst = 0.5*(Ktmp + Ktmp');
        % prevent loss of PSD/symmetry due to numerical errors
end