% SCRIPTPROCESSSABES Example preprocessing data script for [1] for
% cross-validation. 
%
%  [1] ] O Doherty J. et al, "Nonhuman Primate Reaching with Multichannel
%  Sensorimotor Cortex Electrophysiology.” Zenodo, May 26, 2020.
%  doi:10.5281/zenodo.3854034. https://doi.org/10.5281/zenodo.788569
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com


clc;

dirData = '.\..\experimental';
nameData = 'indy_20160927_04.mat';
pathData = fullfile(dirData,nameData);

nFolds = 5;
nameTrain = 'train.mat';
nameTest = 'test.mat';
pathTrain = fullfile(dirData,nameTrain);
pathTest = fullfile(dirData,nameTest);

if ~isfile(pathData)
    error('Check save location of experimental data');
end

if ~isfile(pathTrain)
    fprintf('Generating train and test data\n');
    dataS = load(pathData);
    tlen = length(dataS.t);
    endInds = round(linspace(0,tlen,nFolds+1));
    
    trnInds = 1:endInds(end-1);
    tstInds = endInds(end-1)+1:tlen;
    
    allFields = fieldnames(dataS);
    nFields = length(allFields);
    signalFieldInds = false(1,nFields);
    for fInd = 1:nFields
        field = allFields{fInd};
        if size(dataS.(field),2) == tlen
            signalFieldInds(fInd) = true;
        end
    end
    signalFields = allFields(signalFieldInds);
    copyFields = allFields(~signalFieldInds);
    
    trnS = struct;
    tstS = struct;
    for fInd = 1:length(copyFields)
        field = copyFields{fInd};
        trnS.(field) = dataS.(field);
        tstS.(field) = dataS.(field);
    end
    
    for fInd = 1:length(signalFields)
        field = signalFields{fInd};
        trnS.(field) = dataS.(field)(:,trnInds);
        tstS.(field) = dataS.(field)(:,tstInds);
    end

    save(pathTrain,'-struct','trnS');
    save(pathTest,'-struct','tstS');
else
    fprintf('Found train and test data\n');
end