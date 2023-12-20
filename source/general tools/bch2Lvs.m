function leaves = bch2Lvs(lf,nBch,nLvs)
% Get indices of leaf lf out of nLvs from each of nBch branches
% lf: leaf index off each branch
% nBch: number of branches
% nLvs: number of leaves per branch
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com
leaves = lf + [0:nBch-1]*nLvs;
end

