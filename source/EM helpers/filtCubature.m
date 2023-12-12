function [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtCubature(x,K,n,wts,dsig,A,Q,as,bs)
% FILTCUBATURE Decode latent states using cubature filter
%  One step of the cubature filter of [1].
%  INPUTS:
%  x: previous updated latent state
%  K: covariance of previous estimate
%  n: current Poisson observation
%  wts: cubature weights (generated in getCubePts.m)
%  dsig: cubatutre points (generated in getCubePts.m)
%  A: dynamics matrix
%  Q: latent noise covariance
%  as: baseline firing alphas
%  bs: modulation depth betas
%
%  [1] Song C.Y., Shanechi M.M., "Unsupervised learning of stationary and 
%      switching dynamical sytsem models from Poisson observations",  
%      Journal of Neural Engineering, Dec. 2023
%
%  Author: Christian Song, June 2023, song.christian.y(at)gmail(dot)com

    xPrd = A*x;
    Ktmp = A*K*A' + Q;
    KPrd = 0.5*(Ktmp + Ktmp');
    KPrdi = KPrd^-1;
    
    if ~any(isfinite(n))
        xUpd = xPrd;
        KUpd = KPrd;
    else
        L = chol(KPrd,'lower');
        sigPts = L*dsig + xPrd;
        if any(isfinite(n))
            cifdels = exp(as' + bs'*sigPts);  %dimNt x nPts     
            muNcX = cifdels; %dimNt x nPts 
            covNcX = muNcX; %dimNt x nPts 
            
            muNcN = muNcX*wts'; % E[nt|n1:t-1], dimNt x 1
            mumuNcN = (wts.*muNcX)*muNcX'; % dimNt x dimNt
            Ktmp = diag(covNcX*wts') + mumuNcN - muNcN*muNcN';
            KNN = 0.5*(Ktmp + Ktmp');
            KNNi = KNN^-1;
            KXN = (wts.*(sigPts-xPrd)) * (muNcX - muNcN)';

            xUpd = xPrd + KXN*KNNi*(n-muNcN);
            Ktmp = KPrd - KXN*KNNi*KXN';
            KUpd = 0.5*(Ktmp+Ktmp');
        end
    end

    try
        chol(KUpd,'lower');
    catch
        [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtLocalLap(A,Q,as,bs,x,K,n);
    end
end

