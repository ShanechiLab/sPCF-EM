function bailCheck = estepSEM(theta,nt,prmS,memO,sumO,funO)
% ESTEPSEM E-step of switch EM, smooth brain and regime states
%  Estimate brain and regime states for a swiching dynamical system for
%  switch EM parameterized by theta from Poisson 
%  observations. Requires 3 objects initialized in fitSwitchEM.
%  See [1,2] for details and derivations.
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
%     smoothFunc:   str for choosing which switching smoother to use.
%                   Choose between 'song'(default),'kim'. See [1] for
%                   details.
%  memO: MemSEM object for preallocated memory for intermediate products
%  sumO: SumSEM object for collecting and summing calculated expected value 
%        terms for M-step.
%  funO: FuncOrg object for stationary filter selection, decoded regime
%        probability calculation, and misc functions
%  OUTPUT:
%  bailCheck: bool for if unexpected error interrupted estimation. No other
%             expected values are explicitly returned and instead are
%             collected in sumO. See smoothStates for standalone smoothing
%             function.
% 
%  [1] Song C. et al, "Modeling and Inference Methods for Switching
%      Regime-Dependent Dynamical Systems with Multiscale Neural
%      Observations",  Journal of Neural Engineering, Oct. 2022
%  [1] Song C.Y., Shanechi M.M., "Unsupervised learning of stationary and 
%      switching dynamical sytsem models from Poisson observations",  
%      Journal of Neural Engineering, Dec. 2023
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    if ~exist('prmS','var')
        prmS = struct;
    end

    fields = {'smoothFunc'};
    defaults = {'song'};
    prmS = insertStructDefaults(prmS,fields,defaults);

    [~,tlen] = size(nt);

    dimXt = length(theta.x0Mean);
    dimSt = length(theta.sInit);
    
    smtPars = {'song','kim'};
    smtList = {@smtM,@smtB};
    smtSel = find(contains(smtPars,lower(prmS.smoothFunc)),1);
    smtHand = smtList{smtSel};
    
    decHand = @decM;

    %% Decoding
    bailCheck = false;
    try
        decBase();
        for t = 2:tlen
            KDec = decHand(t);
            if ~isfinitereal(KDec)
                bailCheck = true;
                break;
            end
        end
    catch
        bailCheck = true;
    end
    memO.xSmtBch(:,:) = memO.xMats(:,tlen,:,1);
    memO.KSmtBch(:,:,:) = memO.KMats(:,:,tlen,:,1);
    PSmt = memO.PStDec(:,tlen);
    
    if bailCheck == false
        for t = tlen-1:-1:1
            smtHand(t);
            if ~isfinitereal(memO.KSmtBch)
                bailCheck = true;
                break;
            end
        end
        smtBase();
    end
    
    %% Nested Functions
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
            
            memO.xMats(:,tTo,j,1) = xUpd;
            memO.xMats(:,tTo,j,2) = theta.x0Mean;
            memO.KMats(:,:,tTo,j,1) = KUpd;
            memO.KMats(:,:,tTo,j,2) = theta.x0Cov;
            memO.KMats(:,:,tTo,j,3) = KPrdi;

            % memO.PNtCondSt(:,j) = getExpctdCif(xPrd,KPrd,j);
            memO.PNtCondSt(:,j) = 1 - funO.noSpikeProb(xPrd,KPrd,j);
        end
        PStUpdPre = exp(centWithUpper(logPStUpd,100));
        
        if isfinitereal(PStUpdPre) && dimSt > 1
            PStUpd = (1/sum(PStUpdPre))*PStUpdPre;
        else
            PStUpd = theta.sInit;
        end
        memO.PStDec(:,tTo) = PStUpd;
        
        [mu,~] = mixtureSS(reshape(...
                       memO.xMats(:,tTo,:,1),dimXt,[]),...
                       reshape(memO.KMats(:,:,tTo,:,1),...
                       dimXt,dimXt,[]),...
                       memO.PStDec(:,tTo));
                   
        %% Prediction
        if all(isfinite(nt(:,tTo)))
            PNtPrd = memO.PNtCondSt*theta.sInit;
        else
            PNtPrd = [];
        end

        %% Packing
        sumO.addToDec(mu,PStUpd,PNtPrd);
    end % 0 -> 1

    function sigma = decM(tTo) % t-1 -> t
        tF = tTo-1;
        
        %% Mixing
        PStMixPre = theta.sTran' .* memO.PStDec(:,tF);
        sumT = 1./sum(PStMixPre,1);
        if ~all(isfinitereal(sumT))
            PStMixPre = exp(centWithUpper(log(PStMixPre),1));
            sumT = 1./sum(PStMixPre,1);
        end
        PStMix = sumT .* PStMixPre; 
            % column j: P(S_t-1 | S_t=j,H_t-1)
        
        for j = 1:dimSt
            [mu,sigma] = mixtureSS(reshape(...
                           memO.xMats(:,tF,:,1),dimXt,[]),...
                           reshape(memO.KMats(:,:,tF,:,1),...
                           dimXt,dimXt,[]),...
                           PStMix(:,j));
            memO.xTmpBch(:,j) = mu; % x_t-1|t-1 (S_t=m)
            memO.KTmpBch(:,:,j) = sigma;
        end
        
        %% Stepping
        logPStUpd = zeros(dimSt,1);
        PStPrd = theta.sTran * memO.PStDec(:,tF);
        for j = 1:dimSt
            [xUpd,KUpd,xPrd,KPrd,KPrdi] = funO.applyFilt(j,...
                                    memO.xTmpBch(:,j),...
                                    memO.KTmpBch(:,:,j),...
                                    nt(:,tTo));
            if dimSt > 1
                logPStUpd(j) = real(funO.calcLogPStUpd(...
                                  PStPrd(j),j,xPrd,xUpd,...
                                  KPrdi,KUpd,...
                                  nt(:,tTo)));
            else
                logPStUpd(j) = 0;
            end
            
            memO.xMats(:,tTo,j,1) = xUpd;
            memO.xMats(:,tTo,j,2) = memO.xTmpBch(:,j);
            memO.KMats(:,:,tTo,j,1) = KUpd;
            memO.KMats(:,:,tTo,j,2) = memO.KTmpBch(:,:,j);
            memO.KMats(:,:,tTo,j,3) = KPrdi;

            % memO.PNtCondSt(:,j) = getExpctdCif(xPrd,KPrd,j);
            memO.PNtCondSt(:,j) = 1 - funO.noSpikeProb(xPrd,KPrd,j);
        end
        PStUpdPre = exp(centWithUpper(logPStUpd,100));
        
        %% Estimation
        if isfinitereal(PStUpdPre) && dimSt > 1
            PStUpd = (1/sum(PStUpdPre))*PStUpdPre;
        else
            PStUpd = PStPrd;
        end
        memO.PStDec(:,tTo) = PStUpd;
        [mu,sigma] = mixtureSS(reshape(...
                       memO.xMats(:,tTo,:,1),dimXt,[]),...
                       reshape(memO.KMats(:,:,tTo,:,1),...
                       dimXt,dimXt,[]),...
                       memO.PStDec(:,tTo));
        
        %% Prediction
        if all(isfinite(nt(:,tTo)))
            PNtPrd = memO.PNtCondSt*PStPrd; 
                %P(N_t >= 1| H_t-1)
        else
            PNtPrd = [];
        end
        
        %% Packing
        sumO.addToDec(mu,PStUpd,PNtPrd);
    end

    function smtM(tTo)
        tF = tTo + 1;
        for j = 1:dimSt
            A = theta.Acell{j}; Q = theta.Qcell{j};
            xPrdMix = A*memO.xMats(:,tF,j,2);
            KPrdMix = A*memO.KMats(:,:,tF,j,2)*A'+Q;
            
            memO.xTmpBch(:,j) = memO.xSmtBch(:,j) - xPrdMix;
            memO.KTmpBch(:,:,j) = memO.KSmtBch(:,:,j) - KPrdMix;
        end
        
        %% Direct EM Terms
        sumO.incField('PStSmtSum',PSmt); % PSmt from tF
        for j = 1:dimSt
            Jright = theta.Acell{j}'*memO.KMats(:,:,tF,j,3);
            [xSmtEM,KSmtEM,Jj] = SmoothStepGen(memO.xMats(:,tF,j,2),...
                                               memO.KMats(:,:,tF,j,2),...
                                               memO.xTmpBch(:,j),...
                                               memO.KTmpBch(:,:,j),...
                                               Jright);
            KDifEM = memO.KSmtBch(:,:,j)*Jj';
            
            sumO.incXYFields(j,memO.xSmtBch(:,j),xSmtEM,...
                             memO.KSmtBch(:,:,j),KDifEM,KSmtEM,...
                             PSmt(j));
            if isfinite(nt(1,tF))
                cifOver = funO.checkOver(j,memO.xSmtBch(:,j),...
                                         memO.KSmtBch(:,:,j));
                sumO.addToTot(memO.xSmtBch(:,j),...
                              memO.KSmtBch(:,:,j),PSmt(j),...
                              nt(:,tF),cifOver);
            end
        end
        
        %% Regime Smoothing
        PStMixPre = theta.sTran' .* memO.PStDec(:,tTo);
        PStMix = (1./sum(PStMixPre,1)) .* PStMixPre;
            % column j: P(S_t-1 | S_t=j,H_t-1)
        jointPre = PStMix .* PSmt';
        PStJoint = (1/sum(jointPre(:))) .* jointPre;
        PSmt = sum(PStJoint,2);
        
        sumO.incField('sTranNumSum',PStJoint);
        sumO.incField('sTranDenSum',PSmt);
        
        %% Branching
        for j = 1:dimSt
            Jright = theta.Acell{j}'*memO.KMats(:,:,tF,j,3);
            for i = 1:dimSt
                lfInd = indBch(j,dimSt,i);
                [xEst,KEst] = SmoothStepGen(memO.xMats(:,tTo,i,1),...
                                            memO.KMats(:,:,tTo,i,1),...
                                            memO.xTmpBch(:,j),...
                                            memO.KTmpBch(:,:,j),...
                                            Jright);
                memO.xTmpBch2(:,lfInd) = xEst;
                memO.KTmpBch2(:,:,lfInd) = KEst;
            end
        end
        
        %% Condensation
        for i = 1:dimSt
            if (1/PSmt(i)) ~= Inf
                jointPStQuot = (1/PSmt(i)) * PStJoint(i,:)';
            else
                jointPStQuot = PSmt;
            end
            leafInds = bch2Lvs(i,dimSt,dimSt);
            [mu,sigma] = mixtureSS(reshape(...
                           memO.xTmpBch2(:,leafInds),dimXt,[]),...
                           reshape(memO.KTmpBch2(:,:,leafInds),...
                           dimXt,dimXt,[]),...
                           jointPStQuot);
            memO.xSmtBch(:,i) = mu;
            memO.KSmtBch(:,:,i) = sigma;
        end
    end


    function smtB(tTo)

        tF = tTo + 1;
        %% Direct EM Terms
        sumO.incField('PStSmtSum',PSmt); % PSmt from tF
        for j = 1:dimSt
            A = theta.Acell{j}; Q = theta.Qcell{j};
            xDecTo = memO.xMats(:,tF,j,2);
            KDecTO = memO.KMats(:,:,tF,j,2);
            [xSmtEM,KSmtEM,Jj] = SmoothStepGen(xDecTo,KDecTO,...
                                               [],[],[],A,Q,...
                                               memO.xSmtBch(:,j),...
                                               memO.KSmtBch(:,:,j));
            KDifEM = memO.KSmtBch(:,:,j)*Jj';
            sumO.incXYFields(j,memO.xSmtBch(:,j),xSmtEM,...
                             memO.KSmtBch(:,:,j),KDifEM,KSmtEM,...
                             PSmt(j));
            if isfinite(nt(1,tF))
                cifOver = funO.checkOver(j,memO.xSmtBch(:,j),...
                                         memO.KSmtBch(:,:,j));
                sumO.addToTot(memO.xSmtBch(:,j),...
                              memO.KSmtBch(:,:,j),PSmt(j),...
                              nt(:,tF),cifOver);
            end
        end
        
        %% Regime Smoothing
        PStMixPre = theta.sTran' .* memO.PStDec(:,tTo);
        PStMix = (1./sum(PStMixPre,1)) .* PStMixPre;
            % column j: P(S_t-1 | S_t=j,H_t-1)
        jointPre = PStMix .* PSmt';
        PStJoint = (1/sum(jointPre(:))) .* jointPre;
        PSmt = sum(PStJoint,2);
        
        sumO.incField('sTranNumSum',PStJoint);
        sumO.incField('sTranDenSum',PSmt);
        
        
        %% Branching
        for j = 1:dimSt
            A = theta.Acell{j}; Q = theta.Qcell{j};
            
            for i = 1:dimSt
                [xEst,KEst] = SmoothStepGen(memO.xMats(:,tTo,i,1),...
                                            memO.KMats(:,:,tTo,i,1),...
                                            [],[],[],A,Q,...
                                            memO.xSmtBch(:,j),...
                                            memO.KSmtBch(:,:,j));
                lfInd = indBch(j,dimSt,i);
                memO.xTmpBch2(:,lfInd) = xEst;
                memO.KTmpBch2(:,:,lfInd) = KEst;
            end
        end
        
        %% Condensation
        for i = 1:dimSt
            if (1/PSmt(i)) ~= Inf
                jointPStQuot = (1/PSmt(i)) * PStJoint(i,:)';
            else
                jointPStQuot = PSmt;
            end
            leafInds = bch2Lvs(i,dimSt,dimSt);
            [mu,sigma] = mixtureSS(reshape(...
                           memO.xTmpBch2(:,leafInds),dimXt,[]),...
                           reshape(memO.KTmpBch2(:,:,leafInds),...
                           dimXt,dimXt,[]),...
                           jointPStQuot);
            memO.xSmtBch(:,i) = mu;
            memO.KSmtBch(:,:,i) = sigma;
        end
        
    end

    function smtBase
        tF = 1;
        sumO.incField('PStSmtSum',PSmt);
        for j = 1:dimSt
            A = theta.Acell{j}; Q = theta.Qcell{j};
            [xEst,KEst,Jj] = SmoothStepGen(theta.x0Mean,theta.x0Cov,...
                                          [],[],[],A,Q,...
                                          memO.xSmtBch(:,j),...
                                          memO.KSmtBch(:,:,j));
            memO.xTmpBch(:,j) = xEst;
            memO.KTmpBch(:,:,j) = KEst;
            KDifEM = memO.KSmtBch(:,:,j)*Jj';
            
            sumO.incXYFields(j,memO.xSmtBch(:,j),...
                             memO.xTmpBch(:,j),...
                             memO.KSmtBch(:,:,j),KDifEM,...
                             memO.KTmpBch(:,:,j),...
                             PSmt(j));
            if isfinite(nt(1,tF))
                cifOver = funO.checkOver(j,memO.xSmtBch(:,j),...
                                         memO.KSmtBch(:,:,j));
                sumO.addToTot(memO.xSmtBch(:,j),...
                              memO.KSmtBch(:,:,j),PSmt(j),...
                              nt(:,tF),cifOver);
            end
        end
        
        [x0Smt,K0Smt] = mixtureSS(reshape(...
                                  memO.xTmpBch,dimXt,[]),...
                                  reshape(memO.KTmpBch,...
                                  dimXt,dimXt,[]),...
                                  PSmt);
        sumO.incField('x0SmtSum',x0Smt);
        sumO.incField('K0SmtSum',K0Smt);
        sumO.incField('sInitSum',PSmt);
    end

    function expctdCif = getExpctdCif(xPrd,KPrd,s)
        a = theta.alphaCell{s};
        b = theta.betaCell{s};
        bKb = 0.5 * sum(b' .* (b' * KPrd),2);
        
        expctdCif = exp(a' + b'*xPrd + bKb);
    end
end