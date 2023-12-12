function prmS = insertStructDefaults(prmS,fields,defaults,filename)
    if isempty(prmS)
        prmS = struct;
    end
    if ~exist('filename','var')
        filename = [];
    end
    if length(filename) <= 0
        print = false;
    else
        print = true;
    end
    for fInd = 1:length(fields)
        field = fields{fInd};
        if ~isfield(prmS,field)
            prmS.(field) = defaults{fInd};
            if print
                fprintf('%s: Using default value for %s\n',filename,field);
            end
        end
    end
end