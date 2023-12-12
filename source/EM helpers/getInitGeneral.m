function thetaInit = getInitGeneral(trnS,prmS)
% GETINITGENERAL Parameter initializaion for EM
%  General function for getting an intial set of model parameters (theta) 
%  for a stationary/switching dynamical system. Will try to load a provided
%  theta, otherwise will randomly generate using genInitThetaData.
%  INPUTS:
%  trnS: struct containing training data and pertinent info
%     FIELDS: see fitSwitchEM
%  prmS: struct containing initialization / loading settings
%     FIELDS (main):
%     initGiven:    str for how to try loading initial theta. Choose
%                   between [] (default) for no loading, 'prms' for
%                   checking in prmS, or 'saved' for looking in
%                   prmS.initLocation
%     initLocation: path for saved initial theta if initGiven set to
%                   'saved'
%  OUTPUT:
%  thetaNew: struct containing initial system parameters
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    fields = {'initGiven','initLocation'};
    default = {[],[]};
    prmS = insertStructDefaults(prmS,fields,default);
    
    needInit = true;
    %% Try to load in a valid theta
    if ~isempty(prmS.initGiven)
        loadS = [];
        switch lower(prmS.initGiven)
            case 'true'
                if isfield(trnS,'thetaTrue')
                    loadS = trnS.thetaTrue;
                elseif isfield(prmS,'thetaTrue')
                    loadS = prmS.thetaTrue;
                end
            case 'prms'
                if isfield(prmS,'thetaInit')
                    loadS = prmS.thetaInit;
                elseif isfield(prmS,'theta')
                    loadS = prmS.theta;
                end
            case 'saved'
                pathInit = prmS.initLocation;
                if isfile(pathInit)
                    loadS = load(pathInit);
                end
        end
        
        if iscell(loadS)
            loadS = loadS{1};
        end
        if isfield(loadS,'Acell')
            thetaInit = loadS;
            needInit = false;
        elseif isfield(loadS,'thetaInit')
            if isfield(loadS.thetaInit,'Acell')
                thetaInit = loadS.thetaInit;
                needInit = false;
            end
        elseif isfield(loadS,'theta')
            if isfield(loadS.theta,'Acell')
                thetaInit = loadS.theta;
                needInit = false;
            end
        elseif isfield(loadS,'thetaCell')
            if isfield(loadS.thetaCell{1},'Acell')
                thetaInit = loadS.thetaCell{1};
                needInit = false;
            end
        end
    end
    
    %% If failed to load or init randomly
    if needInit
        baseS = prmS;
        baseS.initGiven = [];

        thetaInit = genInitThetaData(baseS,trnS);
    end

end