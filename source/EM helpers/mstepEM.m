function thetaNew = mstepEM(sumO,prmS,thetaOld)
% MSTEPSEM M-step of switch EM, smooth brain and regime states
%  Estimate brain and regime states for a swiching dynamical system for
%  switch EM parameterized by theta from Poisson observations. Requires 3
%  objects initialized in fitSwitchEM. See [1] for details.
%
%  INPUTS:
%  sumO: SumSEM object with collected expected value terms from E-step.
%  prmS: struct of estimation settings
%     FIELDS (main):
%     statObs:  bool for fixed observation parameters across regimes
%     statDyn:  bool for fixed dynamic parameters across regimes
%     mWeight:  double for how much to mix in previous theta into new
%               theta. 1(default) for no mixing, 0 for no new theta
%     mDiagQ:   str for whether to contrain Q matrices to diagonals. Choose
%               between 'none'(default) or 'all'.
%  thetaOld: struct containing system parameters from previous iteration
%  OUTPUT:
%  thetaNew: struct containing new system parameters
% 
%  [1] Song C.Y., Shanechi M.M., "Unsupervised learning of stationary and 
%      switching dynamical sytsem models from Poisson observations",  
%      Journal of Neural Engineering, Dec. 2023
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    if ~exist('prmS','var')
        prmS = struct;
    end
    fields = {'mDiagQ','statObs',...
              'mWeight','statDyn','rescaleInit'};
    defaults = {'none',false,...
                1,false,1};
    prmS = insertStructDefaults(prmS,fields,defaults);
    
    dimXt = size(thetaOld.Acell{1},1);
    dimSt = length(thetaOld.sInit);
    dimNt = size(thetaOld.betaCell{1},2);
    
    xDiffSum = sumO.xDiffSum;
    xAuM1Sum = sumO.xAuM1Sum;
    xAutoSum = sumO.xAutoSum;
    PStSmtSum = sumO.PStSmtSum;
    x0SmtSum = sumO.x0SmtSum;
    K0SmtSum = sumO.K0SmtSum;
    sTranNumSum = sumO.sTranNumSum;
    sTranDenSum = sumO.sTranDenSum;
    sInitSum = sumO.sInitSum;
    nSegs = sumO.nSegs;

    xTot = sumO.xTot;
    KTot = sumO.KTot;
    PTot = sumO.PTot;
    ntTot = sumO.ntTot;

    alphaCur = thetaOld.alphaCell;
    betaCur = thetaOld.betaCell;

    maximizeFunc = @maximizePoiLikeli;
    optimOpt = optimoptions('fminunc',...
                            'Algorithm','trust-region',... 
                            'SpecifyObjectiveGradient',true,...
                            'HessianFcn','objective',... 
                            'Display','off',...
                            'OptimalityTolerance',1e-5,...
                            'StepTolerance',1e-5,...
                            'MaxIterations',200);


    Acell = cell(dimSt,1); Qcell = cell(dimSt,1); 
    alphaCell = cell(dimSt,1); betaCell = cell(dimSt,1);

    if prmS.statDyn
        A = sum(xDiffSum,3)*(sum(xAuM1Sum,3)^-1);
        Acell(:) = {A};
        
        Qnum = sum(xAutoSum,3) - A*sum(xDiffSum,3)' ...
               - sum(xDiffSum,3)*A' + A*sum(xAuM1Sum,3)*A';
        Qden = sum(PStSmtSum);
        Qtemp = Qnum * Qden^-1;
        switch prmS.mDiagQ
            case 'none'
                Q = 0.5*(Qtemp + Qtemp');
            otherwise
                Q = diag(diag(Qtemp));
        end
        Qcell(:) = {Q};
    else
        for k = 1:dimSt
            A = xDiffSum(:,:,k) * xAuM1Sum(:,:,k)^-1;
            Acell{k} = A;

            Qnum = xAutoSum(:,:,k) - A*xDiffSum(:,:,k)' ...
                   - xDiffSum(:,:,k)*A' + A*xAuM1Sum(:,:,k)*A';
            Qden = PStSmtSum(k);
            Qtemp = Qnum * Qden^-1;
            switch prmS.mDiagQ
                case 'all'
                    Qcell{k} = diag(diag(Qtemp));
                otherwise
                    Qcell{k} = 0.5*(Qtemp + Qtemp');
            end
        end
    end

    if prmS.statObs == false
        for k = 1:dimSt
            inds = sumO.incInds & (sumO.regInds==k);
            [aNew,bNew] = maximizeFunc(ntTot(:,inds),...
                                xTot(:,inds),PTot(:,inds),...
                                KTot(:,:,inds),...
                                alphaCur{k},betaCur{k},...
                                optimOpt);
            alphaCell{k} = aNew;
            betaCell{k} = bNew;
        end
    else
        inds = sumO.incInds;
        [aNew,bNew] = maximizeFunc(ntTot(:,inds),...
                                xTot(:,inds),PTot(:,inds),...
                                KTot(:,:,inds),...
                                alphaCur{1},betaCur{1},...
                                optimOpt);
        alphaCell(:) = {aNew};
        betaCell(:) = {bNew};
    end

    
    sTranPre = sTranNumSum .*(1./sTranDenSum');
    sInitPre = (1/nSegs).*sInitSum;
    x0CovPre = (1/nSegs).*K0SmtSum;

    thetaNew = struct;
    thetaNew.Acell = Acell;
    thetaNew.Qcell = Qcell;
    thetaNew.alphaCell = alphaCell;
    thetaNew.betaCell = betaCell;
    thetaNew.sTran = (1./sum(sTranPre,1)).* sTranPre;
    
    if prmS.rescaleInit ~= 1
        tmp = (1./sum(sInitPre)) .* sInitPre;
        factor = prmS.rescaleInit;
        thetaNew.sInit = rescaleDist(tmp,factor);
    else
        thetaNew.sInit = (1./sum(sInitPre)) .* sInitPre;
    end
    
    thetaNew.x0Mean = (1/nSegs).*x0SmtSum;
    thetaNew.x0Cov = 0.5*(x0CovPre + x0CovPre');
    

    if isfield(prmS,'mWeight') && exist('thetaOld','var')
        thetaNew = mixTheta(thetaOld,thetaNew,prmS.mWeight);
    end
end