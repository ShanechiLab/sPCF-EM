
function [diff,diffStruct] = getThetaDiff(theta1,theta2)
    field1 = fieldnames(theta1);
    field2 = fieldnames(theta2);
    ignoreFields = {};
    fields = setdiff(intersect(field1,field2),ignoreFields);
    
    
    nFields = length(fields);
    diffStruct = struct;
    diff = 0;
    for fInd = 1:nFields
        field = fields{fInd};
        if iscell(theta1.(field))
            nCell = length(theta1.(field));
            if isnumeric(theta1.(field){1})
                diffStruct.(field) = 0;
                for k = 1:nCell
                    tempDiff = theta1.(field){k}-theta2.(field){k};
                    absDiff = sum(abs(tempDiff(:)));
                    diff = diff + absDiff;
                    diffStruct.(field) = diffStruct.(field) + ...
                                         absDiff;
                end
            end
        elseif isnumeric(theta1.(field))
            diffStruct.(field) = 0;
            tempDiff = theta1.(field)-theta2.(field);
            absDiff = sum(abs(tempDiff(:)));
            diff = diff + absDiff;
            diffStruct.(field) = diffStruct.(field) + absDiff;
        end
    end
end