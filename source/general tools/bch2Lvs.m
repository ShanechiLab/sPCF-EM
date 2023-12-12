function leaves = bch2Lvs(lf,nBch,nLvs)
    % Get indices of leaf lf out of nLvs from each of nBch branches
    % lf: leaf index off each branch
    % nBch: number of branches
    % nLvs: number of leaves per branch
leaves = lf + [0:nBch-1]*nLvs;
end

