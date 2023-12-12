function out = indBch(branch,leaves,leaf)
% indBch: index of leaf across all leaves
% branch: which branch
% leaves: leaves per branch
% leaf: index of leaf on branch
    out = (branch-1)*leaves + leaf;
end

