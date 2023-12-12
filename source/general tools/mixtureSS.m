function [mu,sigma] = mixtureSS(mus,sigmas,weights)
% MIXTURESS Get mean and covariance of mixture of Gaussians
%  mus: dxM means of M d-dimensional Gaussians
%  sigmas: dxdxM covariances of Gaussians
%  weights: M weights of Gaussians
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    if numel(weights) > 1
        [N,M] = size(mus);

        if size(weights,1) == 1
            weights = weights';
        end
        if length(weights) ~= M
            error('weights must be same length as number of Gaussians');
        end

        mu = mus * weights;
        muD = mus - mu;
        
        tmp = reshape(reshape(sigmas,N^2,M)*weights,N,N) ...
                + (weights'.*muD)*muD';
        sigma = 0.5*(tmp+tmp');
    else
        mu = mus;
        sigma = sigmas;
    end
end

