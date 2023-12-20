classdef MemSEM < handle
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    properties
        theta
        
        xMats
        KMats
        PStDec
        xTmpBch
        KTmpBch
        xTmpBch2
        KTmpBch2
        PNtCondSt
        xSmtBch
        KSmtBch
    end
    
    methods
        function obj = MemSEM(theta,dimXt,dimSt,dimNt,...
                              maxLen,prmS)
            tlen = maxLen;
            obj.xMats = zeros(dimXt,tlen,dimSt,2);
            obj.KMats = zeros(dimXt,dimXt,tlen,dimSt,3);
            obj.PStDec = zeros(dimSt,tlen);
            obj.xTmpBch = zeros(dimXt,dimSt);
            obj.KTmpBch = zeros(dimXt,dimXt,dimSt);
            obj.xTmpBch2 = zeros(dimXt,dimSt^2);
            obj.KTmpBch2 = zeros(dimXt,dimXt,dimSt^2);
            obj.PNtCondSt = zeros(dimNt,dimSt);
            obj.xSmtBch = zeros(dimXt,dimSt);
            obj.KSmtBch = zeros(dimXt,dimXt,dimSt);
            
            if ~exist('prmS','var') || isempty(prmS)
                prmS = struct;
            end
            
            fields = {};
            defaults = {};
            prmS = insertStructDefaults(prmS,fields,defaults);
            
            obj.updateTheta(theta);
        end
        
        function obj = updateTheta(obj,theta)
            obj.theta = theta;
        end
    end 
end