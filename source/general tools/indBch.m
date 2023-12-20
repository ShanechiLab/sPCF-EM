function out = indBch(branch,leaves,leaf)
% indBch: index of leaf across all leaves
% branch: which branch
% leaves: leaves per branch
% leaf: index of leaf on branch
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
    out = (branch-1)*leaves + leaf;
end

