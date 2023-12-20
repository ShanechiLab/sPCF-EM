function lrnS = getDefaultParams(dimXtEst,dimStEst)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    lrnS = struct;
    switch nargin
        case 1
            lrnS.dimXtEst = dimXtEst;
        case 2
            lrnS.dimXtEst = dimXtEst;
            lrnS.dimStEst = dimStEst;
    end

    lrnS.statDyn = false; 
    lrnS.statObs = false; 
    
    lrnS.fsYt = 100;

    % additional settings
    % initialization settings for genInitThetaData
    lrnS.initGiven = []; lrnS.initLocation = [];
    lrnS.pStay = 0.9999; lrnS.AdiagVal = 0.99;
    lrnS.initDiagQ = false; lrnS.Qparams = [0.025, 0.002];
    lrnS.initDiagR = false;
    % additional M-step settings
    lrnS.mWeight = 1; lrnS.mDiagQ = 'none'; % for mstepEM
    lrnS.mDiagR = 'none';
    % E-step settings
    lrnS.filtSelect = 'cubeInf'; lrnS.nPts = 5; % for FuncOrg
    lrnS.smoothFunc = 'song'; % for estepSEM
    
    % settings for fitSwitchEM
    lrnS.obsPred = true;
    lrnS.verbose = true;
    lrnS.saveFreq = 3;

end

