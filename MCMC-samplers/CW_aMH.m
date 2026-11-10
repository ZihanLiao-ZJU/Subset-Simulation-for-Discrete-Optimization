classdef CW_aMH
    % Component-wise modified MH in whitened Gaussian-reference coordinates.
    properties (Constant)
        SedSav = false
    end
    properties
        BayStc
    end
    methods
        function [u,y,Ncal,baystc] = SamGen(obj,u_sed,y_sed,Ns)
            baystc = obj.BayStc;
            [mu,L] = GaussianPrior(baystc);
            n = baystc.Ndim; Nc = size(u_sed,2);
            eta = L\(u_sed-mu);
            sigma = max(std(eta,0,2),0.1);
            u = zeros(n,Nc,Ns); y = zeros(size(y_sed,1),Nc,Ns);
            Ncal = 0;
            for k = 1:Ns
                if k==1 && obj.SedSav
                    u(:,:,k)=u_sed; y(:,:,k)=y_sed; continue
                end
                candidate = eta + sigma.*(2*rand(n,Nc)-1);
                reject = log(rand(n,Nc)) > min(0,-0.5*(candidate.^2-eta.^2));
                candidate(reject) = eta(reject);
                proposed = mu+L*candidate;
                moved = any(candidate~=eta,1);
                u(:,:,k)=u_sed; y(:,:,k)=y_sed;
                if any(moved)
                    idx = find(moved);
                    baystc = baystc.EvlY(proposed(:,idx));
                    yc = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                    Ncal = Ncal+baystc.Ncal;
                    ratio = yc(1,:)-y_sed(1,idx);
                    accept = ~isnan(ratio) & log(rand(1,numel(idx)))<=min(0,ratio);
                    take = idx(accept);
                    u(:,take,k)=proposed(:,take); y(:,take,k)=yc(:,accept);
                    eta(:,take)=candidate(:,take);
                end
                u_sed=u(:,:,k); y_sed=y(:,:,k);
            end
            baystc.X=reshape(u,n,[]);
            baystc.Y=reshape(y(3:end,:,:),size(y,1)-2,[]);
            baystc.Ncal=Ncal;
        end
    end
end
