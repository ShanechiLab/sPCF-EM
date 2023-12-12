function [mngr,bailCheck] = filterStates(theta,nt,prmS)
% FILTERSTATES Causally estimate brain and regime states
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
%     filtSelect:   str for choosing stationary filter to embed in switch
%                   filter, can choose 'cubeInf','cube','laplace'
%     saveInter:    bool for whether to collect additional matrices such
%                   that smoothStates can be run afterwards by setting mngr
%                   as a field in the prmS for 'smoothStates.m'
%  OUTPUT:
%  mngr: struct containing filtered state and estimation covariances
%     FIELDS (main):
%     xDec:         dxT decoded x1:T
%     KDec:         dxdxT covariance of decoded xt
%     PStEst:       MxD probability of decoded regimes
%     xPrdOne:      one step ahead prediction of xt
%     yPrdN:        N step ahead prediction of yt where N is period between
%                   subsequent yt
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
    fields = {'saveInter','filtSelect'};
    defaults = {false,'laplace'};
    prmS = insertStructDefaults(prmS,fields,defaults);
    
    [dimNt,tlen] = size(nt);
    x0Mean = theta.x0Mean;
    x0Cov = theta.x0Cov;
    dimXt = length(x0Mean);
    dimSt = length(theta.sInit);
    funO = FuncOrg(theta,prmS);

    decHand = @decM;
    
    saveFlag = prmS.saveInter;
    
    if saveFlag
        xMats = zeros(dimXt,tlen,dimSt,2);       % dec
        KMats = zeros(dimXt,dimXt,tlen,dimSt,3); % dec
    else
        xMats = zeros(dimXt,tlen,dimSt);       % dec
        KMats = zeros(dimXt,dimXt,tlen,dimSt); % dec
    end

    PStEst = zeros(dimSt,tlen); %dec   
    xDec = zeros(dimXt,tlen);
    KDec = zeros(dimXt,dimXt,tlen);
    
    xTmpBch = zeros(dimXt,dimSt);
    KTmpBch = zeros(dimXt,dimXt,dimSt);
    
    xPrdOne = zeros(dimXt,tlen);
    PNtPrd = zeros(dimNt,tlen); % true P(n_t|n_1:t-1)
    xPrdOneCondSt = zeros(dimXt,dimSt);
    PNtCondSt = zeros(dimNt,dimSt);
    
    x0 = zeros(dimXt,1);
    K0 = zeros(dimXt,dimXt);
    x0(:,1) = x0Mean;
    K0(:,:) = x0Cov;
    
    %% Decoding 
    bailCheck = false;
    try
        decBase();
        for t = 2:tlen
            Kcur = decHand(t);
            if ~isfinitereal(Kcur)
                bailCheck = true;
                break;
            end
        end
    catch
        bailCheck = true;
    end
    
    
    %% Clean up
    mngr = struct;
    mngr.xMats = xMats;
    mngr.KMats = KMats;
    mngr.PStEst = PStEst;
    mngr.x0 = x0;
    mngr.K0 = K0;
    mngr.xDec = xDec;
    mngr.KDec = KDec;
    mngr.xPrdOne = xPrdOne;
    mngr.PNtPrd = PNtPrd; % P(n_t|n_1:t-1)
    
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
            if saveFlag
                KMats(:,:,tTo,j,3) = KPrdi;
            end
            
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
        PStEst(:,tTo) = PStUpd;
        
        [mu,sigma] = mixtureSS(reshape(...
                       xMats(:,tTo,:,1),dimXt,[]),...
                       reshape(KMats(:,:,tTo,:,1),...
                       dimXt,dimXt,[]),...
                       PStEst(:,tTo));
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
        PStMixPre = theta.sTran' .* PStEst(:,tF);
        PStMix = (1./sum(PStMixPre,1)) .* PStMixPre; 
            % column j: P(S_t-1 | S_t=j,H_t-1)
        PStPrd = theta.sTran*PStEst(:,tF);
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
            KMats(:,:,tTo,j,1) = KUpd;
            if saveFlag
                xMats(:,tF,j,2) = xTmpBch(:,j);
                KMats(:,:,tF,j,2) = KTmpBch(:,:,j);
                KMats(:,:,tTo,j,3) = KPrdi;
            end
            
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
        PStEst(:,tTo) = PStUpd;
        [mu,sigma] = mixtureSS(reshape(...
                       xMats(:,tTo,:,1),dimXt,[]),...
                       reshape(KMats(:,:,tTo,:,1),...
                       dimXt,dimXt,[]),...
                       PStEst(:,tTo));
        xDec(:,tTo) = mu;
        KDec(:,:,tTo) = sigma;
        
        %% Prediction
        xPrdOne(:,tTo) = xPrdOneCondSt*PStPrd;
            % x_t|t-1 = sum x_t|t-1(S_t) P(S_t|H_t-1)
        if all(isfinite(nt(:,tTo)))
            PNtPrd(:,tTo) = PNtCondSt * PStPrd; %P(N_t | H_t-1)
        end
    end

    function expctdCif = getExpctdCif(xPrd,KPrd,s)
        a = theta.alphaCell{s};
        b = theta.betaCell{s};
        bKb = 0.5 * sum(b' .* (b' * KPrd),2);
        
        expctdCif = exp(a' + b'*xPrd + bKb);
    end
end