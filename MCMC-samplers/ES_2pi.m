classdef ES_2pi < MCMCBase
    % Gaussian-reference elliptical slice sampling, one angle per chain.
    % A shared angle is essential when the Gaussian covariance is correlated.
    properties (Constant)
        Ntry = 1000
    end
    methods
        function obj = ES_2pi(model)
            if nargin < 1, model = []; end
            obj@MCMCBase(model);
        end

        function [u,y,Ncal,baystc] = SamGen(obj,varargin)
            [baystc,u_sed,y_sed,Ncal,saveSeed,L,Nc,Ns] = obj.prepare(varargin{:});
            n = baystc.Ndim;
            u = zeros(n,Nc,Ns);
            y = zeros(baystc.Nfun,Nc,Ns);
            for k=1:Ns
                if k==1 && saveSeed
                    u(:,:,k)=u_sed; y(:,:,k)=y_sed; continue
                end
                nu=L*randn(n,Nc);
                threshold=y_sed(1,:)+log(rand(1,Nc));
                angle=2*pi*rand(1,Nc);
                lo=angle-2*pi; hi=angle;
                active=1:Nc;
                for attempt=1:obj.Ntry
                    proposed=u_sed(:,active).*cos(angle(active))+nu(:,active).*sin(angle(active));
                    baystc=baystc.EvlY(proposed);
                    yc=[baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                    Ncal=Ncal+baystc.Ncal;
                    accept=isfinite(yc(1,:)) & yc(1,:)>=threshold(active);
                    accepted=active(accept);
                    u(:,accepted,k)=proposed(:,accept); y(:,accepted,k)=yc(:,accept);
                    active=active(~accept);
                    if isempty(active), break; end
                    positive=angle(active)>0;
                    hi(active(positive))=angle(active(positive));
                    lo(active(~positive))=angle(active(~positive));
                    angle(active)=lo(active)+rand(1,numel(active)).*(hi(active)-lo(active));
                end
                if ~isempty(active)
                    error('SeMIS:SliceBracket','Elliptical slice bracket did not find an accepted point.');
                end
                u_sed=u(:,:,k); y_sed=y(:,:,k);
            end
            baystc.X=reshape(u,n,[]);
            baystc.Y=reshape(y(3:end,:,:),size(y,1)-2,[]);
            baystc.Ncal=Ncal;
        end
    end
end
