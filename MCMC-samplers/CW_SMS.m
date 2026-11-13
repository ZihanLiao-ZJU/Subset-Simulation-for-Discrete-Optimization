classdef CW_SMS < MCMCBase
    % Component-wise ensemble stretch MH. Updates complementary chain groups
    % sequentially, so each proposal partner is fixed during its MH update.
    properties (Constant)
        A = 3
    end
    methods
        function obj = CW_SMS(model)
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
                if Nc==1
                    % A one-chain stretch ensemble is undefined. Symmetric
                    % Gaussian random walk is a valid one-chain MH fallback.
                    proposal=u_sed+L*randn(n,1);
                    baystc=baystc.EvlY(proposal);
                    yc=[baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                    Ncal=Ncal+baystc.Ncal;
                    ratio=sum(yc(1:2,:))-sum(y_sed(1:2,:));
                    if ~isnan(ratio) && log(rand)<=min(0,ratio)
                        u_sed=proposal; y_sed=yc;
                    end
                else
                    order=randperm(Nc);
                    groups={order(1:ceil(Nc/2)),order(ceil(Nc/2)+1:end)};
                    for group=1:2
                        active=groups{group}; partners=groups{3-group};
                        for dim=1:n
                            partner=partners(randi(numel(partners),1,numel(active)));
                            z=((sqrt(obj.A)-1/sqrt(obj.A))*rand(1,numel(active))+1/sqrt(obj.A)).^2;
                            proposal=u_sed(:,active);
                            proposal(dim,:)=u_sed(dim,partner)+z.*(u_sed(dim,active)-u_sed(dim,partner));
                            baystc=baystc.EvlY(proposal);
                            yc=[baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                            Ncal=Ncal+baystc.Ncal;
                            % One-dimensional stretch has Jacobian z^(1-1).
                            ratio=sum(yc(1:2,:),1)-sum(y_sed(1:2,active),1);
                            accept=~isnan(ratio) & log(rand(1,numel(active)))<=min(0,ratio);
                            u_sed(:,active(accept))=proposal(:,accept);
                            y_sed(:,active(accept))=yc(:,accept);
                        end
                    end
                end
                u(:,:,k)=u_sed; y(:,:,k)=y_sed;
            end
            baystc.X=reshape(u,n,[]);
            baystc.Y=reshape(y(3:end,:,:),size(y,1)-2,[]);
            baystc.Ncal=Ncal;
        end
    end
end
