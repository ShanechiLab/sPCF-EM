function thetaInit = genInitThetaData(prmS,trnS)
% GENINITTHETADATA Parameter initializaion for EM based on data
%  Randomly generates initial system parameters based on provided data and
%  initialization settings.
%  INPUTS:
%  prmS: struct containing initialization / loading settings
%     FIELDS (required):
%     dimXtEst:     d, dimension of latent brain state of model
%     dimStEst:     M, # of regimes of model
%     pStay:        diagonal values of regime transition matrix
%     sTran:        prespecified regime transition matrix (alternative to
%                   using dimStEst and pStay)
%     sTranEst:     same as sTran
%     FIELDS (optional):
%     Qparams:      1x2 double of mean and variance of Gaussian that Q's
%                   eigenvalues are drawn from
%     initDiagQ:    bool for initializing Q as diagonal matrices
%     statObs:      bool for fixed observation parameters across regimes
%     AdiagVal:     double for what to initialize diagonal values of A as,
%                   defaulted to 1 (identity matrices)
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    if isfield(prmS,'dimXtEst')
        dimXt = prmS.dimXtEst;
    else
        error('Need dimXtEst field in prmS');
    end
    fs = trnS.fs;
    if isfield(prmS,'sTran')
        dimSt = size(prmS.sTran,1);
        sTran = prmS.sTran;
    elseif isfield(prmS,'dimStEst') && isfield(prmS,'pStay')
        dimSt = prmS.dimStEst;
        sTran = generateMarkovParams(prmS.pStay,dimSt);
    elseif isfield(prmS,'sTranEst')
        dimSt = size(prmS.sTranEst,1);
        sTran = prmS.sTranEst;
    elseif isfield(prmS,'dimStEst')
        dimSt = prmS.dimStEst;
        pStay = 0.9999;
        sTran = generateMarkovParams(pStay,dimSt);
    else
        error('Need dimStEst (number of regimes) in prmS');
    end
    
    fields = {'Qparams','initDiagQ',...
			  'statObs','AdiagVal'};
	defaults = {[0.025,0.002],false,...
				false,0.9};
    prmS = insertStructDefaults(prmS,fields,defaults);
    
    %% Prepping data

    obsNtAll = prepInputObs(trnS,prmS);
    ntTot = catCell(obsNtAll);
    tTot = size(ntTot,2);

    dimNt = size(obsNtAll{1},1);
    ntSum = sum(ntTot,2);
    
    thetaInit.Acell = cell(dimSt,1);
    thetaInit.Qcell = cell(dimSt,1);
    thetaInit.alphaCell = cell(dimSt,1);
    thetaInit.betaCell = cell(dimSt,1);
    thetaInit.x0Mean = zeros(dimXt,1);
    thetaInit.x0Cov = eye(dimXt);
    thetaInit.sInit = (1/dimSt).*ones(dimSt,1);
    thetaInit.sTran = (1./sum(sTran,1)).*sTran;

    %% Generating A Matrices
    thetaInit.Acell(:) = {prmS.AdiagVal.*eye(dimXt)};

    
    %% Generating xt Noise Matrix
    Qprms = prmS.Qparams;
    for m = 1:dimSt
        thetaInit.Qcell{m} = getRandCov(dimXt,Qprms(1),Qprms(2),...
                                        'norm',prmS.initDiagQ);
    end
    QvarSum = zeros(dimXt,1);
    for m = 1:dimSt
        QvarSum = QvarSum + diag(thetaInit.Qcell{m});
    end
    Qvar = 1/dimSt .* QvarSum;
    
    %% Generating Poisson Parameters
    baseFires = ntSum ./ (tTot/fs); %avg firerate
    
    baseL = 0.8.*(0.5.*baseFires);
    baseH = 1.2.*(0.5.*baseFires);
    thetaInit.alphaCell(:) = {zeros(1,dimNt)};
    for c = 1:dimNt
        bases = getUnifRand(baseL(c),baseH(c),1,dimSt);
        for m = 1:dimSt
            a = log(bases(m)) - log(fs);
            thetaInit.alphaCell{m}(c) = a;
        end
    end
    
    xMax = 4*sqrt(1000*max(Qvar)); % V[xt] ~ 1/(1-decay)^2 * V[wt]
    betaRngs = abs(log(60./min(baseFires,59)))./(dimXt*xMax);
    
    if prmS.statObs
        alphaTemp = thetaInit.alphaCell{1};
        thetaInit.alphaCell(:) = {alphaTemp};
        
        betaTemp = zeros(dimXt,dimNt);
        for c = 1:dimNt
            bL = -1*betaRngs(c);
            bH = betaRngs(c);
            betaTemp(:,c) = getUnifRand(bL,bH,dimXt,1);
        end
        thetaInit.betaCell(:) = {betaTemp};
    else
        for m = 1:dimSt
            betaTemp = zeros(dimXt,dimNt);
            for c = 1:dimNt
                bL = -1*betaRngs(c);
                bH = betaRngs(c);
                betaTemp(:,c) = getUnifRand(bL,bH,dimXt,1);
            end
            thetaInit.betaCell{m} = betaTemp;
        end
    end

end