function thetaMix = mixTheta(thetaOld,thetaNew,pt)

    thetaMix = thetaNew;
    if pt < 1
        mixFields_cell = {'Acell'};
        mixFields_value = {'sTran','sInit','x0Cov','x0Mean'};
        dimSt = length(thetaNew.Acell);
        for fInd = 1:length(mixFields_cell)
            field = mixFields_cell{fInd};
            if isfield(thetaNew,field)
                for k = 1:dimSt
                    thetaMix.(field){k} = (1-pt)*thetaOld.(field){k} ...
                                            + pt*thetaNew.(field){k};
                end
            end
        end

        for fInd = 1:length(mixFields_value)
            field = mixFields_value{fInd};
            if isfield(thetaNew,field)
                thetaMix.(field) = (1-pt)*thetaOld.(field) ...
                                     + pt*thetaNew.(field);
            end
        end
    end
end