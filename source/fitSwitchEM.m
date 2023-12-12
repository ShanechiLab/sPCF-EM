function mdlS = fitSwitchEM(trnS,lrnPrmS,nIter,pathModel)
% FITSWITCHEM Fit system parameters using switch EM
%  Learn parameters for a switching dynamical system with Poisson
%  observations using a switch EM framework. See [1] for details and
%  derivations.
% 
%  INPUTS:
%  trnS: struct containing training data and pertinent info
%     FIELDS (required):
%     fs:       sampling frequency (Hz), double
%     obsNt:    Poisson observations, formatted as single matrix
%               N(# of channels) x T(# of samples) or cell array of
%               separate trials with consistent N but allows for varying T.
%               obsNt not required if Gaussian observations provided
%     FIELDS (optional):
%     stateXt:  behavior observations, formatted similar to obsNt as 
%               d(# of dimensions) x either T or T+1. T+1 for when x0
%               provided. Does not impact learned parameters. Only used to
%               track decoding performance of learned parameters through a
%               similarity transform from decoded xt to stateXt.
%     stateSt:  regime observations formatted the same as obsNt with
%               columns representing either 1D true regime labels within
%               [1,M] or equivalent MD one hot represenations. Does not
%               impact learned parameters. Only used to track decoding
%               performance of learned parameters.
%  lrnPrmS: struct containing learn settings
%     FIELDS (main):
%     dimXtEst: d, dimension of latent brain state of model
%     dimStEst: M, # of regimes of model
%     obsPred:  bool for whether to track observation prediction metrics
%     statObs:  bool for fixed observation parameters across regimes
%     statDyn:  bool for fixed dynamic parameters across regimes
%     verbose:  bool for printing tracked metrics
%     filtSelect: str for choosing stationary filter to embed in switch
%               filter, can choose 'cubeInf','cube','laplace', see
%               EM helpers/FuncOrg.m'
%     Additional fields: see 'EM helpers/getDefaultParams.m'
%  nIter: # of EM iterations to run
%  pathModel: path to where learned model is saved. If provided, is also 
%             where progress is saved to stop/resume learning. Increasing
%             the # of iterations can also pick up from last saved
%             iteration
% 
%  OUTPUT:
%  mdlS: struct containing learned model parameters and any tracked metrics
%     FIELDS (guaranteed):
%     thetaCell:cell array 1 x nIter+1 of learned parameters from each EM
%               iteration with thetaCell{1} being the initial parameters
%               and thetaCell{end} being the final learned parameters
%     static:   struct containing learning parameters and misc info if EM
%               iterations unexpectedly fail
%     FIELDS (metrics):
%     ntPrdMet: PP = 2*AUC-1, One step ahead prediction of spikes nt using
%               P(nt|nt-1)
%     xtCC:     Correlation coefficient between provided stateXt and
%               decoded xt post similarity transform
%     stAcc:    Accuracy of decoded st against provided stateSt
% 
%  [1] Song C.Y., Shanechi M.M., "Unsupervised learning of stationary and 
%      switching dynamical sytsem models from Poisson observations",  
%      Journal of Neural Engineering, Dec. 2023
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    %% Setting setup
    xtGiven = isfield(trnS,'stateXt');
    stGiven = isfield(trnS,'stateSt');
    if ~exist('lrnPrmS','var') || isempty(lrnPrmS)
        lrnPrmS = struct;
    end
    if ~isfield(lrnPrmS,'verbose')
        verbose = false;
        curFile = '';
    else
        verbose = lrnPrmS.verbose;
        curFile = 'fitSwitchEM';
    end
    varList = {'obsPred','verbose','saveFreq'};
    defaultVals = {true,true,3};
    lrnPrmS = insertStructDefaults(lrnPrmS,varList,defaultVals,curFile);

    if isfield(lrnPrmS,'printName')
        msg = sprintf(', %s',lrnPrmS.printName);
    else
        msg = '';
    end
    if verbose
        fprintf('%s%s',curFile,msg);
    end
    
    %% Metrics to track
    genFields = {'runTime'};
    genFSizes = 3;
    obsFields = {'ntPrdMet'};
    obsFSizes = 1;
    xFields = {'xtCC'};
    xFSizes = 1;
    sFields = {'stAcc'};
    sFSizes = 1;
    
    mdlFields = genFields;
    mdlFSizes = genFSizes;
    
    if lrnPrmS.obsPred
        mdlFields = [mdlFields, obsFields];
        mdlFSizes = [mdlFSizes, obsFSizes];
    end
    if xtGiven
        mdlFields = [mdlFields, xFields];
        mdlFSizes = [mdlFSizes, xFSizes];
    end
    if stGiven
        mdlFields = [mdlFields, sFields];
        mdlFSizes = [mdlFSizes, sFSizes];
    end

    %% Checking for and Loading Existing Model
    if exist('pathModel','var') && ~isempty(pathModel)
        savePerIter = true;
        dirLearn = fileparts(pathModel);
        if ~exist(dirLearn,'dir')
            mkdir(dirLearn);
        end
    else
        savePerIter = false;
    end
    if savePerIter && exist(pathModel,'file')
        [mdlS,newMdl,toRun,startIter,newIter] = loadModel(pathModel,...
                                                  lrnPrmS,nIter);
        if newIter ~= nIter
            nIter = newIter;
        end
    else
        newMdl = true;
        toRun = true;
        startIter = 1;
    end

    if toRun == false
        if verbose
            fprintf(' Existing completed detected,');
        end
    else
        %% Start prep for EM             
        if verbose
            fprintf('Start: %s',datetime);
        end
        
        if newMdl == true
            %% Initialization
            thetaInit = getInitGeneral(trnS,lrnPrmS);

            mdlS = struct;
            mdlS.thetaCell = cell(1,nIter+1);
            mdlS.static = struct;
            for fInd = 1:length(mdlFields)
                field = mdlFields{fInd};
                mdlS.(field) = zeros(mdlFSizes(fInd),nIter);
            end
            mdlS.thetaCell{1} = thetaInit;
            mdlS.static.lrnPrmS = lrnPrmS;
        end
        %% Prepping observations and objects
        obsNtAll = prepInputObs(trnS,lrnPrmS);
        nSegs = length(obsNtAll);
        if lrnPrmS.obsPred
            ntTot = catCell(obsNtAll);
        end
        if xtGiven == true
            xtWhole = prepXtCat(trnS);
        else
            xtWhole = [];
        end
        if stGiven == true
            stWhole = catCell(trnS.stateSt);
        end
        [sumClct,memMngr,funMngr] = initObjects(obsNtAll,...
                                                mdlS.thetaCell{startIter},...
                                                lrnPrmS.obsPred,...
                                                xtGiven,stGiven,lrnPrmS);

        %% EM
        bailedOut = false;
        iter = startIter;
        if verbose
            fprintf(', iter:');
        end
        while iter <= nIter && ~bailedOut
            if verbose
                fprintf('%d,',iter);
            end
            tic;
            thetaEst = mdlS.thetaCell{iter};
            
            funMngr.resetTheta(thetaEst);
            memMngr.updateTheta(thetaEst);
            sumClct.reset();
            %% E-Step
            for s = 1:nSegs
                obsNtSeg = obsNtAll{s};
                segBail = estepSEM(thetaEst,obsNtSeg,lrnPrmS,...
                                      memMngr,sumClct,funMngr);
                bailedOut = bailedOut || segBail;
            end
            mdlS.runTime(1,iter) = toc; tic;

            if bailedOut
                mdlS.static.bailedOut = true;
                mdlS.static.bailedOutInd = iter;
                mdlS.static.bailedOutTheta = thetaEst;
                iter = iter - 1;
            else
                %% M-Step
                thetaNew = mstepEM(sumClct,lrnPrmS,thetaEst);
                mdlS.runTime(2,iter) = toc; tic;
                
                %% metrics
                if verbose
                    fprintf('(');
                end
                if lrnPrmS.obsPred
                    spPP = getPP(sumClct.PNtPrdTot,ntTot);
                    mdlS.ntPrdMet(:,iter) = spPP;
                    if verbose
                        fprintf('n%0.3g',spPP);
                    end
                end
                if xtGiven == true
                    simTran = learnProjSignal(xtWhole,sumClct.xDecTot);
                    mdlS.thetaCell{iter}.simTran = simTran;
                    xCC = averageCorrCoef(xtWhole,simTran*sumClct.xDecTot);
                    mdlS.xtCC(iter) = xCC;
                    if verbose
                        fprintf('x%0.3g',xCC);
                    end
                end
                if stGiven == true
                    stAcc = max(calcStAccuracy(sumClct.PStDecTot,stWhole));
                    mdlS.stAcc(iter) = stAcc;
                    if verbose
                        fprintf('s%0.3g',stAcc);
                    end
                end
                if verbose
                    fprintf('pSty%0.3g', mean(diag(thetaEst.sTran)));
                end

                mdlS.runTime(3,iter) = toc; tic;
                mdlS.thetaCell{iter + 1} = thetaNew;
                if verbose
                    fprintf('),');
                end
            end
            iter = iter + 1;
            if savePerIter
                saveCond = lrnPrmS.saveFreq <= 1 ...
                           || mod(iter-1,lrnPrmS.saveFreq)==1 ...
                           || ~(iter <= nIter);
                if saveCond
                    save(pathModel,'-struct','mdlS','-v7.3');
                end
            end
        end
        
        if bailedOut
            mdlS = trimModel(mdlS,iter - 1);
            if savePerIter
                save(pathModel,'-struct','mdlS','-v7.3');
            end
            if verbose
                fprintf(', bailed out. ');
            end
        end
    end
    if verbose
        fprintf(' End: %s \n', datetime);
    end
end


function [sumO,memO,funO] = initObjects(obsCell,theta,...
                                        obsPred,xtGiven,stGiven,lrnS)
    dimXt = size(theta.Acell{1},1);
    dimSt = length(theta.sInit);
    dimNt = length(theta.alphaCell{1});
    nSegs = length(obsCell);
    segLens = zeros(1,nSegs);
    for s = 1:nSegs
        segLens(s) = size(obsCell{s},2);
    end
    sviseFlag = false;
    funO = FuncOrg(theta,lrnS);
    sumO = SumSEM(dimXt,dimSt,dimNt,segLens,...
                  obsPred,xtGiven,stGiven,sviseFlag);
    memO = MemSEM(theta,dimXt,dimSt,dimNt,max(segLens),lrnS);
end
