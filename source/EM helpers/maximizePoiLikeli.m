function [alphaNew,betaNew] = maximizePoiLikeli(nt,xEst,PStEst,KEst,alphaCur,betaCur,optimOpt)
% MAXIMIZEPOILIKLI M-step update for Poisson parameters
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    [dimNt,tlen] = size(nt);
    dimXt = size(xEst,1);
    
    if isempty(PStEst)
        PStEst = ones(1,tlen);
    end
    
    if ~exist('optimOpt','var')
        optimOpt = optimoptions('fminunc',...
                            'Algorithm','trust-region',... 
                            'SpecifyObjectiveGradient',true,...
                            'Display','off');
    end
    if iscell(optimOpt)
        optimOpt = optimOpt{1};
    end
    
    guesses = [alphaCur;betaCur];
    alphaNew = zeros(1,dimNt);
    betaNew = zeros(dimXt,dimNt);

    if isempty(KEst)
        for c = 1:dimNt
            zEst = fminSmpX(guesses(:,c),optimOpt,...
                            nt(c,:),xEst,PStEst);
            alphaNew(:,c) = zEst(1);
            betaNew(:,c) = zEst(2:end);
        end
    else
        for c = 1:dimNt
            zEst = fminExpX(guesses(:,c),optimOpt,...
                            nt(c,:),xEst,KEst,PStEst);
            alphaNew(:,c) = zEst(1);
            betaNew(:,c) = zEst(2:end);
        end
    end
end

function zEst = fminExpX(z0,option,nt,xSmt,KSmt,PStSmt)
    dimXt = size(xSmt,1);
    pstnt = PStSmt.*nt;
    try
        zEst = fminunc(@minFunc,z0,option);
        if ~(all(isfinite(zEst)) && all(isreal(zEst)))
            error('');
        end
    catch
        zEst = z0;
    end

    function [g,gradG,hessG] = minFunc(z)
        if size(z,1) ==1
            z = z';
        end
        a = z(1);
        b = z(2:end);
        Kb = reshape(b' * reshape(KSmt,dimXt,[]),dimXt,[]);

        bx = b' * xSmt;
        pexpabx = PStSmt.*exp(a + bx + 0.5*b'*Kb);

        g = sum(pexpabx - pstnt.*(a + bx));

        gradG = zeros(length(z),1);
        gradG(1) =  sum(pexpabx - pstnt);
        gradG(2:end) = sum(pexpabx.*(xSmt + Kb) ...
                                      - pstnt.*xSmt,2);
        xKb = xSmt + Kb;
        htemp = KSmt + reshape(xKb,dimXt,1,[])...
                      .*reshape(xKb,1,dimXt,[]);
        hb = sum(htemp.*reshape(pexpabx,1,1,[]),3);

        ha = sum(pexpabx);
        hab = sum(pexpabx.*(xSmt + Kb),2)';
        hessG = zeros(length(z),length(z));
        hessG(1,2:end) = hab;
        hessG(2:end,1) = hab';
        hessG(1,1) = ha;
        hessG(2:end,2:end) = hb;
    end
end

function zEst = fminSmpX(z0,option,nt,xSmt,PStSmt)
    pstnt = PStSmt.*nt;
    try
        zEst = fminunc(@minFunc,z0,option);
        if ~(all(isfinite(zEst)) && all(isreal(zEst)))
            error('');
        end
    catch
        zEst = z0;
    end

    function [g,gradG,hessG] = minFunc(z)
        if size(z,1) ==1
            z = z';
        end
        a = z(1);
        b = z(2:end);

        bx = b' * xSmt;
        pexpabx = PStSmt.*exp(a + bx);

        g = sum(pexpabx - pstnt.*(a + bx));

        gradG = zeros(length(z),1);
        gradG(1) =  sum(pexpabx - pstnt);
        gradG(2:end) = sum((pexpabx - pstnt).*xSmt,2);

        pExpX = pexpabx.*xSmt;
        hb = pExpX*xSmt';
        ha = sum(pexpabx);        
        hab = sum(pExpX,2)';

        hessG = zeros(length(z),length(z));
        hessG(1,2:end) = hab;
        hessG(2:end,1) = hab';
        hessG(1,1) = ha;
        hessG(2:end,2:end) = hb;
    end
end