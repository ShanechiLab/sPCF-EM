function [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtCubeInf(x,K,n,wts,dsig,A,Q,as,bs)
% FILTCUBEINF Decode latent states using cubature filter (information)
%  One step of the information form of the cubature filter of [1].
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
        Kesum = zeros(size(Q));
        xsum = zeros(size(x));
        
        if all(isfinite(n))
            L = chol(KPrd,'lower');
            sigPts = L*dsig + xPrd;
            cifdels = exp(as' + bs'*sigPts);  %dimNt x nPts     
            muNcX = cifdels; %E[nt|xt=sigpts], dimNt x nPts 
            covNcX = cifdels; %V[nt|xt=sigpts], dimNt x nPts 
            
            muNcN = muNcX*wts'; % E[nt|n1:t-1], dimNt x 1
            
            KXN = (wts.*(sigPts-xPrd)) * (muNcX - muNcN)';
            Ct = KPrdi*KXN;
            nRi = diag(1./(covNcX*wts')); %E[V[nt|xt]|n1:t-1] inv
            G = Ct*nRi;
%             Ktmp = KPrdi*KXN*nRi*KXN'*KPrdi';
            Ktmp = G*Ct';
            Kesum = Kesum + 0.5*(Ktmp+Ktmp');
            xsum = G*(n-muNcN);
        end
        
        KUpdi = KPrdi + Kesum;
        Ktmp = KUpdi^-1;
        KUpd = 0.5*(Ktmp+Ktmp'); 
            % prevent loss of PSD/symmetry due to numerical errors
        xUpd = xPrd + KUpd*xsum;
        
        try
            chol(KUpd,'lower');
        catch
            [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtLocalLap(A,Q,as,bs,x,K,n);
        end
    end
end