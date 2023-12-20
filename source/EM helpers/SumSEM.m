classdef SumSEM < handle
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    properties
        xDiffSum
        xAuM1Sum
        xAutoSum
        PStSmtSum
        x0SmtSum
        K0SmtSum
        sTranNumSum
        sTranDenSum
        sInitSum
        nSegs
        xTot% [seg1_reg1 ... seg1_regM | ... | segN_reg1 ... segN_regM]
            % [seg1: t1r1 ... t1rM t2r1 ... t2rM ... | ... | segN: t1r1 ...]
        KTot
        PTot
        ntTot
        incInds % [1 1 1 0 0  1 ... | ... | 1 1 1 1 0 1 1 ...]
        regInds % [1 2 3 ... M 1 2 .. M | ... | 1 ... M 1 ... M ...]
        PNtPrdTot
        xDecTot
        PStDecTot
        tCur
        tmCur

        predNt = false
        xtGiven = false
        stGiven = false
        tlen
    end
    
    methods
        function obj = SumSEM(dimXt,dimSt,dimNt,...
                              segLens,...
                              obsPred,xtGiven,stGiven,...
                              sviseFlag,prmS)
            if ~exist('sviseFlag','var') || isempty(sviseFlag)
                sviseFlag = false;
            end
            if ~exist('prmS','var')
                prmS = struct;
            end
            fields = {};
            defaultVals = {};
            prmS = insertStructDefaults(prmS,fields,defaultVals);
            
            obj.xDiffSum = zeros(dimXt,dimXt,dimSt);
            obj.xAuM1Sum = zeros(dimXt,dimXt,dimSt);
            obj.xAutoSum = zeros(dimXt,dimXt,dimSt);
            obj.PStSmtSum = zeros(dimSt,1);
            obj.x0SmtSum = zeros(dimXt,1);
            obj.K0SmtSum = zeros(dimXt,dimXt);
            obj.sTranNumSum = zeros(dimSt,dimSt);
            obj.sTranDenSum = zeros(dimSt,1);
            obj.sInitSum = zeros(dimSt,1);
            obj.nSegs = length(segLens);
            obj.tCur = 1;
            obj.tmCur = 1;
            tlen = sum(segLens);
            obj.tlen = tlen;

            if sviseFlag
                combLen = tlen;
                obj.PTot = ones(1,combLen);
                obj.regInds = zeros(1,tlen);
            else
                combLen = dimSt*tlen;
                obj.PTot = zeros(1,combLen);
                obj.regInds = repmat(1:dimSt,1,tlen);
            end
            
            obj.xTot = zeros(dimXt,combLen);
            obj.KTot = zeros(dimXt,dimXt,combLen);
            
            obj.ntTot = zeros(dimNt,combLen);
            obj.incInds = true(1,combLen);

            if obsPred
                obj.PNtPrdTot = zeros(dimNt,tlen);
                obj.predNt = true;
            end
            if xtGiven
                obj.xDecTot = zeros(dimXt,tlen);
                obj.xtGiven = true;
            end
            if stGiven
                obj.PStDecTot = zeros(dimSt,tlen);
                obj.stGiven = true;
            end
        end

        function addToDec(obj,xDec,PStDec,PNtPrd)
            t = obj.tCur;
            if obj.predNt
                obj.PNtPrdTot(:,t) = PNtPrd;
            end
            if obj.xtGiven
                obj.xDecTot(:,t) = xDec;
            end
            if obj.stGiven
                obj.PStDecTot(:,t) = PStDec;
            end
            obj.tCur = obj.tCur + 1; %= tlen+1 at end
        end
        
        function addToTot(obj,xSmt,KAuto,PStSmt,nt,cifOver)
            t = obj.tmCur;
            obj.xTot(:,t) = xSmt;
            obj.KTot(:,:,t) = KAuto;
            obj.PTot(:,t) = PStSmt;
            obj.ntTot(:,t) = nt;

            if cifOver
                obj.incInds(t) = false;
            end
            obj.tmCur = obj.tmCur + 1;
        end

        function incXYFields(obj,s,xSmt,xExt,KAuto,KDiff,KAuM1,...
                             PStSmt)
            xDiff = PStSmt*(KDiff + xSmt*xExt');
            xAuM1 = PStSmt*(KAuM1 + xExt*xExt');
            xAuto = PStSmt*(KAuto + xSmt*xSmt');
            
            obj.incFieldCondSt('xDiffSum',xDiff,s);
            obj.incFieldCondSt('xAuM1Sum',xAuM1,s);
            obj.incFieldCondSt('xAutoSum',xAuto,s);
        end
        
        function incField(obj,field,value)
            obj.(field) = obj.(field) + value;
        end
        
        function incFieldCondSt(obj,field,value,s)
            obj.(field)(:,:,s) = obj.(field)(:,:,s) + value;
        end

        function reset(obj)
            statProps = {'nSegs','incInds','regInds',...
                         'tCur','tmCur','predNt',...
                         'xtGiven','stGiven','tlen'};

            props = setdiff(properties(obj),statProps);
            for pInd = 1:length(props)
                prop = props{pInd};
                obj.(prop)(:) = 0;
            end
            obj.tCur = 1;
            obj.tmCur = 1;
            obj.incInds(:) = true;
        end
    end
    
end

