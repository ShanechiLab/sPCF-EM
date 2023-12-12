classdef FuncOrg < handle
    properties
        filtSelect
        theta
        
        dimSt
        
        wts = []
        cubeDsig = []

        wtsPP = []
        ptsPP = []
    end
    
    methods
        function obj = FuncOrg(theta,prmS)
            
            if ~exist('prmS','var') || isempty(prmS)
                prmS = struct;
            end
            
            fields = {'filtSelect','nPts'};
            defaults = {'cubeInf',5};
            prmS = insertStructDefaults(prmS,fields,defaults);
            
            dimSt = length(theta.sInit);
            obj.dimSt = dimSt;
            obj.resetTheta(theta);
            
            dimXt = length(theta.Acell{1});

            if contains(prmS.filtSelect,'cube')
                [wts,dsig] = getCubePts(dimXt,prmS.nPts);
                obj.wts = wts;
                obj.cubeDsig = dsig;
            else
                switch prmS.filtSelect
                    case 'laplace'
                        
                    otherwise
                        error('unrecognized filtSelect');
                end
            end
            obj.filtSelect = prmS.filtSelect;
            [obj.wtsPP,obj.ptsPP] = getCubePts(dimXt,5);
        end
        
        function obj = reInit(obj)
            
        end
        
        function obj = resetTheta(obj,theta)
            obj.theta = theta;
        end
        
        function cif = getCif(obj,s,x)
            cif = exp(obj.theta.alphaCell{s}' ...
                      + obj.theta.betaCell{s}'*x);
        end
        function over = checkOver(obj,s,xSmt,KSmt)
            over = false;
            axb = obj.theta.alphaCell{s}' + ...
                  obj.theta.betaCell{s}'*xSmt;
              
            bxb = 0.5* sum(obj.theta.betaCell{s}'.* ...
                           obj.theta.betaCell{s}'*KSmt,2);
            if any(axb + bxb > 100)
                over = true;
            end
        end

        function [xUpd,KUpd,xPrd,KPrd,KPrdi] = applyFilt(obj,s,x,K,n)
            A = obj.theta.Acell{s}; Q = obj.theta.Qcell{s}; 
            a = obj.theta.alphaCell{s}; b = obj.theta.betaCell{s};
            switch obj.filtSelect
                case 'laplace'
                    [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtLocalLap(...
                        A,Q,a,b,x,K,n);
                case 'cubeInf'
                    [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtCubeInf(...
                        x,K,n,obj.wts,obj.cubeDsig,A,Q,a,b);
                case 'cube'
                    [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtCubature(...
                        x,K,n,obj.wts,obj.cubeDsig,A,Q,a,b);
            end
        end
        
        function logPStUpd = calcLogPStUpd(obj,PStPrd,s,xP,xU,KPi,KU,N)
            if PStPrd ~= 0
                logDetPart = 0.5*(logdet(KU) + logdet(KPi));
                logNtPart = 0;
                if all(isfinite(N))
                    cifDel = obj.getCif(s,xU); 
                    logTemp = N'*log(cifDel) - sum(cifDel);
                    if isfinitereal(logTemp)
                        logNtPart = logTemp;
                    end
                end
                logPrePart = -0.5*(xU-xP)'*KPi*(xU-xP);
                logPStUpd = logDetPart ...
                            + logNtPart + logPrePart ...
                            + log(PStPrd);
            else
                logPStUpd = log(PStPrd);
            end
        end

        function PnoSpk = noSpikeProb(obj,xPrd,KPrd,j)
            cSpk = length(obj.theta.alphaCell{1});
            PnoSpk = zeros(cSpk,1);

            as = obj.theta.alphaCell{j};
            bs = obj.theta.betaCell{j};
            try
                L = chol(KPrd,'lower');
                sigPts = L*obj.ptsPP + xPrd;
                expabx = exp(as' + bs'*sigPts); %cSpk x nPts
                expexpabx = exp(-expabx); %cSpk x nPts
                PnoSpk(:) = nansum(expexpabx.*obj.wtsPP,2);
            end
        end
    end
end

