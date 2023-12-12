function [mngr,bailCheck] = smoothStates(theta,nt,prmS)
% SMOOTHSTATES Acausally estimate brain and regime states
%  Estimate brain and regime states for a swiching dynamical system
%  parameterized by theta from Poisson observations. See [1] for details
%  and derivations
%  NOTATION:
%     M: # of regimes
%     d: xt dimension
%     N: nt dimension
%     T: # of time samples
% 
%  INPUTS:
%  theta: struct containing system parameters
%     FIELDS (required):
%     Acell:        1xM cell array of dxd dynamics matrices
%     Qcell:        1xM cell array of dxd latent noise covariance
%     alphaCell:    1xM cell array of 1xN base firing rates
%     betaCell      1xM cell array of dxN firing rate modulation depths
%     x0Mean        dx1 mean of initial x0
%     x0Cov         dxd covariance of initial x0
%     sInit         Mx1 initial regime distribution P(S1)
%     sTran         MxM regime transition matrix, column i: P(st|st-1=i)
%  nt: Poisson observations formatted as NxT matrix. 
%  prmS: struct of estimation settings
%     FIELDS (main):
%     smoothFunc:   str for choosing which switching smoother to use.
%                   Choose between 'song'(default),'kim','EC'. See [1] for
%                   details.
%     useDecMngr:   bool for using output of filterStates.m to just run
%                   backward step. Requires mngr output to be set as field
%                   of prmS and for useDecMngr to be true.
%  OUTPUT:
%  mngr: struct containing filtered state and estimation covariances
%     FIELDS (main):
%     xSmt:         dxT smoothed x1:T
%     KSmt:         dxdxT covariance of smoothed xt
%     xDec:         dxT decoded x1:T
%     KDec:         dxdxT covariance of decoded xt
%     PStEst:       MxD probability of decoded regimes
%     xPrdOne:      one step ahead prediction of xt
%     PNtPrd:       one step ahead probability of P(nt>=1|n1:t-1)
%  bailCheck: bool for if unexpected error interrupted filtering
% 
%  [1] Song C. et al, "Modeling and Inference Methods for Switching
%      Regime-Dependent Dynamical Systems with Multiscale Neural
%      Observations",  Journal of Neural Engineering, Oct. 2022
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    if ~exist('prmS','var')
        prmS = struct;
    end
    fields = {'smoothFunc','useDecMngr'};
    defaults = {'song',false};
    prmS = insertStructDefaults(prmS,fields,defaults);
    
    [dimNt,tlen] = size(nt);
    x0Mean = theta.x0Mean;
    x0Cov = theta.x0Cov;
    dimXt = length(x0Mean);
    dimSt = length(theta.sInit);
    funO = FuncOrg(theta,prmS);

    smtPars = {'song','kim','EC'};
    smtList = {@smtM,@smtB,@smtBrb};
    smtSel = find(contains(smtPars,lower(prmS.smoothFunc)),1);
    smtHand = smtList{smtSel};
    
    decHand = @decM;
    
    useMngr = false;
    mngrConds = isfield(prmS,'mngr') ...
                && isstruct(prmS.mngr) ...
                && isfield(prmS.mngr,'KMats') ...
                && ~isfield(prmS.mngr,'smtSel');
    if prmS.useDecMngr && mngrConds
        if tlen == size(prmS.mngr.KMats,3)
            if size(prmS.mngr.KMats,5) == 3
                useMngr = true;
            end
        end
    end
    
    if ~useMngr
        xDec = zeros(dimXt,tlen);
        KDec = zeros(dimXt,dimXt,tlen);
        
        xPrdOne = zeros(dimXt,tlen);
        PNtPrd = zeros(dimNt,tlen); % true P(n_t|n_1:t-1)
        xPrdOneCondSt = zeros(dimXt,dimSt);
        PNtCondSt = zeros(dimNt,dimSt);
    else
        xDec = prmS.mngr.xDec;
        KDec = prmS.mngr.KDec;
    end
    
    PStEst = zeros(dimSt,tlen,2); %dec   
    xSmt = zeros(dimXt,tlen);
    KSmt = zeros(dimXt,dimXt,tlen);
    xTmpBch = zeros(dimXt,dimSt);
    KTmpBch = zeros(dimXt,dimXt,dimSt);
	xTmpBch2 = zeros(dimXt,dimSt^2);
	KTmpBch2 = zeros(dimXt,dimXt,dimSt^2);

    xMats = zeros(dimXt,tlen,dimSt,2);       % dec,smt/mix
    KMats = zeros(dimXt,dimXt,tlen,dimSt,3); % dec,smt/mix,prdi

    x0 = zeros(dimXt,2);
    K0 = zeros(dimXt,dimXt,2);
    x0(:,1) = x0Mean;
    K0(:,:,1) = x0Cov;
    
    %% Decoding
    bailCheck = false;
    if ~useMngr
        decBase();
        for t = 2:tlen
            KCur = decHand(t);
            if ~isfinitereal(KCur)
                bailCheck = true;
                break;
            end
        end
    else
        PStEst(:,:,1) = prmS.mngr.PStEst;
        xMats(:) = prmS.mngr.xMats(:);
        KMats(:) = prmS.mngr.KMats(:);
    end
    
    if bailCheck == false
        xMats(:,tlen,:,2) = xMats(:,tlen,:,1);
        KMats(:,:,tlen,:,2) = KMats(:,:,tlen,:,1);
        PStEst(:,tlen,2) = PStEst(:,tlen,1);
        xSmt(:,tlen) = xDec(:,tlen);
        KSmt(:,:,tlen) = KDec(:,:,tlen);
        for t = tlen-1:-1:1
            smtHand(t);
            if ~isfinitereal(KSmt(:,:,t))
                bailCheck = true;
                break;
            end
        end
        smtBase();
    end

    
    
    %% Clean up
    mngr = struct;
    
    mngr.xMats = xMats;
    mngr.KMats = KMats;
    mngr.PStEst = PStEst;
    mngr.smtSel = smtSel;
    mngr.x0 = x0;
    mngr.K0 = K0;
    mngr.xDec = xDec;
    mngr.KDec = KDec;
    mngr.xSmt = xSmt;
    mngr.KSmt = KSmt;

    if useMngr
        copyFields = {'xPrdOne','PNtPrd'};
        for fInd = 1:length(copyFields)
            field = copyFields{fInd};
            if isfield(prmS.mngr,field)
                mngr.(field) = prmS.mngr.(field);
            end
        end
    else
        mngr.xPrdOne = xPrdOne;
        mngr.PNtPrd = PNtPrd; % P(n_t|n_1:t-1)
    end
    
	function decBase
        tTo = 1;
        logPStUpd = zeros(dimSt,1);
        for j = 1:dimSt
            [xUpd,KUpd,xPrd,KPrd,KPrdi] = funO.applyFilt(j,...
                                        theta.x0Mean,theta.x0Cov,...
                                        nt(:,tTo));
            
            if dimSt > 1
                logPStUpd(j) = real(funO.calcLogPStUpd(...
                                  theta.sInit(j),...
                                  j,xPrd,xUpd,KPrdi,KUpd,...
                                  nt(:,tTo)));
            else
                logPStUpd(j) = 0;
            end
            xMats(:,tTo,j,1) = xUpd;
            KMats(:,:,tTo,j,1) = KUpd;
            KMats(:,:,tTo,j,3) = KPrdi;
            
            xPrdOneCondSt(:,j) = xPrd;
            % PNtCondSt(:,j) = getExpctdCif(xPrd,KPrd,j);
            PNtCondSt(:,j) = 1 - funO.noSpikeProb(xPrd,KPrd,j);
        end
        PStUpdPre = exp(centWithUpper(logPStUpd,100));
        
        if isfinitereal(PStUpdPre) && dimSt > 1
            PStUpd = (1/sum(PStUpdPre))*PStUpdPre;
        else
            PStUpd = theta.sInit;
        end
        PStEst(:,tTo,1) = PStUpd;
        
        [mu,sigma] = mixtureSS(reshape(...
                       xMats(:,tTo,:,1),dimXt,[]),...
                       reshape(KMats(:,:,tTo,:,1),...
                       dimXt,dimXt,[]),...
                       PStEst(:,tTo,1));
        xDec(:,tTo) = mu;
        KDec(:,:,tTo) = sigma;

        
        %% Extra Stuff for Predictions
        xPrdOne(:,tTo) = xPrdOneCondSt*theta.sInit;
        if all(isfinite(nt(:,tTo)))
            PNtPrd(:,tTo) = PNtCondSt*theta.sInit;
        end
    end % 0 -> 1    
    
    function sigma = decM(tTo)
        tF = tTo-1;
        
        %% Mixing
        PStMixPre = theta.sTran' .* PStEst(:,tF,1);
        PStMix = (1./sum(PStMixPre,1)) .* PStMixPre; 
            % column j: P(S_t-1 | S_t=j,H_t-1)
        PStPrd = theta.sTran*PStEst(:,tF,1);
        for j = 1:dimSt
            [mu,sigma] = mixtureSS(reshape(...
                           xMats(:,tF,:,1),dimXt,[]),...
                           reshape(KMats(:,:,tF,:,1),...
                           dimXt,dimXt,[]),...
                           PStMix(:,j));
            xTmpBch(:,j) = mu; % x_t-1|t-1 (S_t=m)
            KTmpBch(:,:,j) = sigma;
        end 

        %% Stepping
        logPStUpd = zeros(dimSt,1);
        for j = 1:dimSt
            [xUpd,KUpd,xPrd,KPrd,KPrdi] = funO.applyFilt(j,...
                                            xTmpBch(:,j),...
                                            KTmpBch(:,:,j),...
                                            nt(:,tTo));
            if dimSt > 1
                logPStUpd(j) = real(funO.calcLogPStUpd(...
                                  PStPrd(j),j,xPrd,xUpd,...
                                  KPrdi,KUpd,...
                                  nt(:,tTo)));
            else
                logPStUpd(j) = 0;
            end
            xMats(:,tTo,j,1) = xUpd;
            xMats(:,tF,j,2) = xTmpBch(:,j);
            KMats(:,:,tTo,j,1) = KUpd;
            KMats(:,:,tF,j,2) = KTmpBch(:,:,j);
            KMats(:,:,tTo,j,3) = KPrdi;
            
            xPrdOneCondSt(:,j) = xPrd;
            % PNtCondSt(:,j) = getExpctdCif(xPrd,KPrd,j);
            PNtCondSt(:,j) = 1 - funO.noSpikeProb(xPrd,KPrd,j);
        end
        PStUpdPre = exp(centWithUpper(logPStUpd,100));
        
        %% Estimation
        if isfinitereal(PStUpdPre) && dimSt > 1
            PStUpd = (1/sum(PStUpdPre))*PStUpdPre;
        else
            PStUpd = PStPrd;
        end
        PStEst(:,tTo,1) = PStUpd;
        [mu,sigma] = mixtureSS(reshape(...
                       xMats(:,tTo,:,1),dimXt,[]),...
                       reshape(KMats(:,:,tTo,:,1),...
                       dimXt,dimXt,[]),...
                       PStEst(:,tTo,1));
        xDec(:,tTo) = mu;
        KDec(:,:,tTo) = sigma;
        
        %% Predicition
        xPrdOne(:,tTo) = xPrdOneCondSt*PStPrd;% x_t|t-1 = sum x_t|t-1(S_t) P(S_t|H_t-1)
        if all(isfinite(nt(:,tTo)))
            PNtPrd(:,tTo) = PNtCondSt * PStPrd; %P(N_t | H_t-1)
        end
    end

    function smtBase
        tF = 1;
        xDecTo = theta.x0Mean;
        KDecTo = theta.x0Cov;
        for j = 1:dimSt
            A = theta.Acell{j}; Q = theta.Qcell{j};
            [xTmp,Ktmp] = SmoothStepGen(xDecTo,KDecTo,...
                                        [],[],[],A,Q,...
                                        xMats(:,tF,j,2),...
                                        KMats(:,:,tF,j,2));
            xTmpBch(:,j) = xTmp;
            KTmpBch(:,:,j) = 0.5*(Ktmp+Ktmp');
        end
        [mu,sigma] = mixtureSS(xTmpBch,KTmpBch,PStEst(:,1,2));
        x0(:,2) = mu;
       	K0(:,:,2) = sigma;
    end

    function smtM(tTo)
        tF = tTo + 1;
        %% MIX prep
        for j = 1:dimSt
            A = theta.Acell{j}; Q = theta.Qcell{j};
            xPrdMix = A*xMats(:,tTo,j,2);
            KPrdMix = A*KMats(:,:,tTo,j,2)*A'+Q;
            xTmpBch(:,j) = xMats(:,tF,j,2) - xPrdMix;
            KDiff = KMats(:,:,tF,j,2) - KPrdMix;
            KTmpBch(:,:,j) = KDiff;
        end
        
        %% Regime Smoothing
        PStMixPre = theta.sTran' .* PStEst(:,tTo,1);
        PStMix = (1./sum(PStMixPre,1)) .* PStMixPre;
            % column j: P(S_t-1 | S_t=j,H_t-1)
        jointPre = PStMix .* PStEst(:,tF,2)';
        PStJoint = (1/sum(jointPre(:))) .* jointPre;
        PStEst(:,tTo,2) = sum(PStJoint,2);
        
        %% Branching/Smoothing
        for j = 1:dimSt
            Jright = theta.Acell{j}'*KMats(:,:,tF,j,3);
            for i = 1:dimSt
                lfInd = indBch(j,dimSt,i);
                [xEst,KEst] = SmoothStepGen(xMats(:,tTo,i,1),...
                                            KMats(:,:,tTo,i,1),...
                                            xTmpBch(:,j),...
                                            KTmpBch(:,:,j),...
                                            Jright);
                xTmpBch2(:,lfInd) = xEst;
                KTmpBch2(:,:,lfInd) = KEst;
            end
        end
        
        %% Condensation
        for i = 1:dimSt
            dP = 1/PStEst(i,tTo,2);
            if dP ~= Inf
                jointPStQuot = dP * PStJoint(i,:)';
            else
                jointPStQuot = PStEst(:,tF,2);
            end
            leafInds = bch2Lvs(i,dimSt,dimSt);
            [mu,sigma] = mixtureSS(xTmpBch2(:,leafInds),...
                                   KTmpBch2(:,:,leafInds),...
                                   jointPStQuot);
            xMats(:,tTo,i,2) = mu; %overwriting saved mix vals
            KMats(:,:,tTo,i,2) = sigma;
        end
        [mu,sigma] = mixtureSS(reshape(...
                        xMats(:,tTo,:,2),dimXt,[]),...
                        reshape(KMats(:,:,tTo,:,2),...
                        dimXt,dimXt,[]),...
                        PStEst(:,tTo,2));
        xSmt(:,tTo) = mu;
        KSmt(:,:,tTo) = sigma;
    end

    function smtB(tTo)
        tF = tTo + 1;
        
        %% Regime Smoothing
        PStMixPre = theta.sTran' .* PStEst(:,tTo,1);
        PStMix = (1./sum(PStMixPre,1)) .* PStMixPre;
            % column j: P(S_t-1 | S_t=j,H_t-1)
        jointPre = PStMix .* PStEst(:,tF,2)';
        PStJoint = (1/sum(jointPre(:))) .* jointPre;
        PStEst(:,tTo,2) = sum(PStJoint,2);
        
        %% Branching
        for j = 1:dimSt
            A = theta.Acell{j}; Q = theta.Qcell{j};
            for i = 1:dimSt
                lfInd = indBch(j,dimSt,i);
                xDecTo = xMats(:,tTo,i,1);
                KDecTo = KMats(:,:,tTo,i,1);
                [xTmp,Ktmp] = SmoothStepGen(xDecTo,KDecTo,...
                                            [],[],[],A,Q,...
                                            xMats(:,tF,j,2),...
                                            KMats(:,:,tF,j,2));
                xTmpBch2(:,lfInd) = xTmp;
                KTmpBch2(:,:,lfInd) = Ktmp;
            end
        end
        
        %% Condensation
        for i = 1:dimSt
            dP = 1/PStEst(i,tTo,2);
            if dP ~= Inf
                jointPStQuot = dP * PStJoint(i,:)';
            else
                jointPStQuot = PStEst(:,tF,2);
            end
            leafInds = bch2Lvs(i,dimSt,dimSt);
            [mu,sigma] = mixtureSS(xTmpBch2(:,leafInds),...
                                   KTmpBch2(:,:,leafInds),...
                                   jointPStQuot);
            xMats(:,tTo,i,2) = mu;
            KMats(:,:,tTo,i,2) = sigma;
        end
        [mu,sigma] = mixtureSS(reshape(...
                        xMats(:,tTo,:,2),dimXt,[]),...
                        reshape(KMats(:,:,tTo,:,2),...
                        dimXt,dimXt,[]),...
                        PStEst(:,tTo,2));
        xSmt(:,tTo) = mu;
        KSmt(:,:,tTo) = sigma;
    end

    function smtBrb(tTo)
        tF = tTo + 1;
        
        %% Branching
        logEcTerm = zeros(dimSt,dimSt);
        for j = 1:dimSt
            A = theta.Acell{j}; Q = theta.Qcell{j};
            for i = 1:dimSt
                lfInd = indBch(j,dimSt,i);
                xPrv = xMats(:,tTo,i,1);
                KPrv = KMats(:,:,tTo,i,1);
                xPrd = A*xPrv;
                Ktmp = A*KPrv*A' + Q;
                KPrd = 0.5*(Ktmp + Ktmp');
                KPrdi = KPrd^-1;
                
                xDif = xMats(:,tF,j,2) - xPrd;
                KDif = KMats(:,:,tF,j,2) - KPrd;
                Jright = A'*KPrdi;
                [xTmp,Ktmp] = SmoothStepGen(xPrv,KPrv,xDif,KDif,Jright);

                xDifEC = xMats(:,tF,j,2) - xPrd;
                logEcTerm(i,j) = -0.5*(xDifEC'*KPrdi*xDifEC ...
                                       + logdet(KPrd,'chol'));

                xTmpBch2(:,lfInd) = xTmp;
                KTmpBch2(:,:,lfInd) = Ktmp;
            end
        end
        
        %% Regime Smoothing
        PStMixPre = theta.sTran' .* PStEst(:,tTo,1);
        PStMix = (1./sum(PStMixPre,1)) .* PStMixPre;
            % column j: P(S_t-1 | S_t=j,H_t-1)
            
        PStMixSmt = zeros(dimSt,dimSt);
            % column j: P(S_t-1 | S_t=j,H_T)
        for j = 1:dimSt
            pTmp = exp(centWithUpper(logEcTerm(:,j),1)).* PStMix(:,j);
            PStMixSmt(:,j) = (1./sum(pTmp)).*pTmp;
        end
        
        jointPre = PStMixSmt .* PStEst(:,tF,2)';
        PStJoint = (1/sum(jointPre(:))) .* jointPre;
        PStEst(:,tTo,2) = sum(PStJoint,2);
        
        %% Condensation
        for i = 1:dimSt
            dP = 1/PStEst(i,tTo,2);
            if dP ~= Inf
                jointPStQuot = dP * PStJoint(i,:)';
            else
                jointPStQuot = PStEst(:,tF,2);
            end
            leafInds = bch2Lvs(i,dimSt,dimSt);
            [mu,sigma] = mixtureSS(xTmpBch2(:,leafInds),...
                                   KTmpBch2(:,:,leafInds),...
                                   jointPStQuot);
            xMats(:,tTo,i,2) = mu;
            KMats(:,:,tTo,i,2) = sigma;
        end
        [mu,sigma] = mixtureSS(reshape(...
                        xMats(:,tTo,:,2),dimXt,[]),...
                        reshape(KMats(:,:,tTo,:,2),...
                        dimXt,dimXt,[]),...
                        PStEst(:,tTo,2));
        xSmt(:,tTo) = mu;
        KSmt(:,:,tTo) = sigma;
    end


    function expctdCif = getExpctdCif(xPrd,KPrd,s)
        a = theta.alphaCell{s};
        b = theta.betaCell{s};
        bKb = 0.5 * sum(b' .* (b' * KPrd),2);
        
        expctdCif = exp(a' + b'*xPrd + bKb);
    end
end