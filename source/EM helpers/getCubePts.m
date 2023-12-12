function [wts, pts] = getCubePts(dimXt,nDeg)
% GETCUBEPTS Get cubature weights and points
%  Generates cubature weights and points for the 5th and 3rd degree
%  spherical-cubature rules of [1]. Defaults to 5th degree, but can also
%  generate 3rd degree with nDeg=3.
%
%  [1] Jia et al, "High-degree cubature Kalman filter", Automatica, 2013
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    if ~exist('nDeg','var')
        nDeg = 5;
    end

    switch nDeg
        case 5
            
            sjP = zeros(dimXt,dimXt*(dimXt-1)/2);
            sjN = zeros(dimXt,dimXt*(dimXt-1)/2);
            ind = 1;
            for l = 2:dimXt
                for k = 1:l-1
                    ek = zeros(dimXt,1); ek(k) = 1;
                    el = zeros(dimXt,1); el(l) = 1;
                    sjP(:,ind) = sqrt(1/2)*(ek + el);
                    sjN(:,ind) = sqrt(1/2)*(ek - el);
                    ind = ind + 1;
                end
            end
            nSj = size(sjP,2);
            nEj = dimXt;
            
            nPts = 2*(dimXt^2) + 1;
            
            wts = zeros(1,nPts);
            pts = zeros(dimXt,nPts);
            
            wts(1) = 2/(dimXt+2);
            pts(:,1) = zeros(dimXt,1);
            ind = 1;
            
            wts(ind+1:ind+nSj) = 1/(dimXt+2)^2;
            pts(:,ind+1:ind+nSj) = sqrt(dimXt+2).*sjP;
            ind = ind + nSj;
            
            wts(ind+1:ind+nSj) = 1/(dimXt+2)^2;
            pts(:,ind+1:ind+nSj) = -sqrt(dimXt+2).*sjP;
            ind = ind + nSj;
            
            wts(ind+1:ind+nSj) = 1/(dimXt+2)^2;
            pts(:,ind+1:ind+nSj) = sqrt(dimXt+2).*sjN;
            ind = ind + nSj;
            
            wts(ind+1:ind+nSj) = 1/(dimXt+2)^2;
            pts(:,ind+1:ind+nSj) = -sqrt(dimXt+2).*sjN;
            ind = ind + nSj;
            
            wts(ind+1:ind+nEj) = (4-dimXt)/(2*(dimXt+2)^2);
            pts(:,ind+1:ind+nEj) = sqrt(dimXt+2).*eye(dimXt);
            ind = ind + nEj;
            
            wts(ind+1:ind+nEj) = (4-dimXt)/(2*(dimXt+2)^2);
            pts(:,ind+1:ind+nEj) = -sqrt(dimXt+2).*eye(dimXt);
        otherwise
            nPts = 2*dimXt;

            wts = (1/nPts).*ones(1,nPts);

            pts = sqrt(dimXt).*[eye(dimXt), -1*eye(dimXt)];
    end
    

    
end