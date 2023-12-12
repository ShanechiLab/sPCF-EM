% SCRIPTLEARNSABES Example learning script with sPCF-EM on real data.
%  Settings and data used to learn the magnitude regime of [1] as
%  visualized in Figure 8. Details on methods can be see in [1]. Data was
%  downloaded from [2] from session '20160927_04'. Velocity of reaches
%  saved to stateXt but only used for learning similarity transforms;
%  sPCF-EM is unsupervised and only uses neural data to learn parameters.
%  Due to runtime considerations, an existing model outputted from method
%  provided.
%
%  [1] Song C.Y., Shanechi M.M., "Unsupervised learning of stationary and 
%  switching dynamical sytsem models from Poisson observations",  
%  Journal of Neural Engineering, Dec. 2023
%
%  [2] ] O Doherty J. et al, "Nonhuman Primate Reaching with Multichannel
%  Sensorimotor Cortex Electrophysiology.” Zenodo, May 26, 2020.
%  doi:10.5281/zenodo.3854034. https://doi.org/10.5281/zenodo.788569
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

clc;
saveModel = true;

%% Loading and path setup
dirData = '.\..\experimental';
nameTrain = 'train.mat';
nameTest = 'test.mat';
pathTrain = fullfile(dirData,nameTrain);
pathTest = fullfile(dirData,nameTest);
pathInit = fullfile(dirData,'thetaInit.mat');
trnS = load(pathTrain);

if saveModel
    pathModel = fullfile(dirData,'model.mat');
else
    pathModel = [];
end

%% Settings
dimXtEst = 10;
dimStEst = 2;
lrnS = getDefaultParams(dimXtEst,dimStEst);
lrnS.statObs = true; 
lrnS.initGiven = 'saved'; lrnS.initLocation = pathInit;

%% Switch EM Learning with PCF
useProvidedModel = false; 
    % learning 300 iterations will take ~2.4 hrs, feel free to use provided
    % model if only interested in method output. 300 was used for [1]. You
    % can also try different # of iterations.  If you're saving the model,
    % future runs will start from the last saved iteration.

nIterations = 300;
if ~useProvidedModel
    mdlEM = fitSwitchEM(trnS,lrnS,nIterations,pathModel);
else
    pathModelProv = fullfile(dirData,'model_provided.mat');
    mdlEM = load(pathModelProv);
end

theta = mdlEM.thetaCell{nIterations};

%% Fitting Similarity Transform / Projection
mngrTrn = filterStates(theta,trnS.obsNt);
xDecTrn = mngrTrn.xDec;
vtTrn = trnS.stateXt;
L = learnProjSignal(vtTrn,[xDecTrn; ones(1,size(xDecTrn,2))]);

%% Test set decoding
tstS = load(pathTest);
nt = tstS.obsNt;
mngr = filterStates(theta,nt);
xDecTst = mngr.xDec;
vDecTst = L*[xDecTst;ones(1,size(xDecTst,2))];
PStDec = mngr.PStEst;

%% Validation
vtTst = tstS.stateXt;
vtCC = averageCorrCoef(vtTst,vDecTst);
fprintf('Decoded Velocity CC: %0.3g\n',vtCC);

%% Regime Plot
tLim = [4977.252,4980.822];
regInd = 1;
cmapRegBlu = [4,79,137]./255;

t = tstS.t;
targOn = t(tstS.targetOnset);
sDec = round(PStDec(regInd,:));

figure;
plot(t,PStDec(regInd,:));
hold on;
for k = 1:length(targOn)
    vertLine(targOn(k),[],'m-');
end
addStPatches(t,sDec,[],cmapRegBlu,0.2);
xlim(tLim);
xlabel('Time(sec)');
ylabel('P(s_t)')
title('Decoded Regime Probability');

