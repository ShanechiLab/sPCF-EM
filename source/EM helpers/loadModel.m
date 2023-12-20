function [mdlS,needNewMdl,needToRun,startIter,newIter] = loadModel(pathModel,prmS,nIter)
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    if ~exist('prmS','var')
        prmS = struct;
    end
    if isfield(prmS,'verbose')
        verbose = prmS.verbose;
    else
        verbose = false;
    end

    mdlS = [];
    needNewMdl = true;
    needToRun = true;
    startIter = 1;
    newIter = nIter;
    if isfile(pathModel)
        loaded = false;
        try
            loadMdl = load(pathModel);
            if isstruct(loadMdl) && isfield(loadMdl,'thetaCell') ...
                                 && isfield(loadMdl,'static')
                loaded = true;
            end
        catch
            if verbose
                fprintf(' failed to load ');
            end
        end
        if loaded
            [done,finIt,newIt] = checkLoadProgress(loadMdl,nIter);
            needNewMdl = false;
            if done
                needToRun = false;
                ldIter = length(loadMdl.thetaCell)-1;
                if finIt < ldIter
                    mdlS = trimModel(loadMdl,finIt);
                    save(pathModel,'-struct','mdlS','-v7.3');
                else
                    mdlS = loadMdl;
                end
            else
                startIter = finIt + 1;
                if ~isempty(newIt)
                    nIter = newIt;
                    mdlS = expandModel(loadMdl,nIter);
                else
                    mdlS = loadMdl;
                end
            end
        end
        if needNewMdl
            [dirPath,filename,ext] = fileparts(pathModel);
            newName = [filename, '_old'];
            backupPath = fullfile(dirPath,[newName, ext]);
            try
                movefile(pathModel,backupPath);
            catch
                if verbose
                    fprintf(' failed to backup ');
                end
            end
        end
    else
        warning('model file doesnt exist');
    end
end

