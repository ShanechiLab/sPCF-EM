function [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtLocalLap(A,Q,as,bs,x,K,N)
% FILTCUBATURE Decode latent states using local Laplace approximation
%  One step of the point process filter [1-2].
%  INPUTS:
%  x: previous updated latent state
%  K: covariance of previous estimate
%  N: current Poisson observation
%  A: dynamics matrix
%  Q: latent noise covariance
%  alpha: baseline firing alphas
%  beta: modulation depth betas
%
%  [1] Eden et al, "Dynamic Analysis of Neural Encoding by Point Process
%      Adaptive Filtering",  Neural Computation, 2004
%  [2] Hsieh et al, "Multiscale modeling and decoding algorithms for
%      spike-field activity", Journal of Neural Engineering, 2019
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    xPrd = A*x;
    Ktmp = A*K*A' + Q;
    KPrd = 0.5*(Ktmp + Ktmp');
    KPrdi = KPrd^-1;
    
    Kesum = zeros(size(Q));
    xsum = zeros(size(x));
    if all(isfinite(N))
        cifdel = exp(as' + bs'*xPrd);
        xtemp =  bs * (N - cifdel);
        Ktemp = (cifdel' .* bs) * bs';
        if isfinitereal(xtemp) && isfinitereal(Ktemp)
            xsum = xtemp;
            Kesum = Ktemp;
        end
    end
    
    KupdTemp = (KPrdi + Kesum)^-1;     
    
    if any(isnan(KupdTemp(:)))
        KUpd = KPrd;
    else
        KUpd = 0.5*(KupdTemp + KupdTemp');
    end
    
    xUpd = xPrd + KUpd*xsum;

end

