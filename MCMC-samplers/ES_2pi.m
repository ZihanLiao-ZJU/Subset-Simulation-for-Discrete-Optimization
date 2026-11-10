classdef ES_2pi
    % Gaussian-reference elliptical slice sampling, one angle per chain.
    % A shared angle is essential when the Gaussian covariance is correlated.
    properties (Constant)
        SedSav = false
        Ntry = 1000
    end
    properties
        BayStc
    end
    methods
        function obj = ConfigureISDPolicy(obj,defaults)
            % SeMCS validates field names and values before this hook. The
            % slice kernel moves every seed and has no aCS scale, carry or
            % acceptance-rate adaptation; those policies are inapplicable.
            if isfield(defaults,'SeedMode') && ~strcmp(defaults.SeedMode,'move')
                error('SeMIS:SamplerDefaults', ...
                    'Elliptical slice sampling supports SeedMode=move only.');
            end
            if isfield(defaults,'MixingMode') && ~strcmp(defaults.MixingMode,'off')
                error('SeMIS:SamplerDefaults', ...
                    'Elliptical slice sampling does not implement aCS mixing diagnostics or score adaptation.');
            end
        end

        function [u,y,Ncal,baystc] = SamGen(obj,u_sed,y_sed,Ns)
            if ~isscalar(Ns)
                error('SeMIS:RaggedSampler', ...
                    'Elliptical slice sampling requires one common chain length.');
            end
            validateattributes(Ns,{'numeric'},{'real','integer','positive','finite'});
            baystc = obj.BayStc;
            [mu,L] = GaussianPrior(baystc);
            n = baystc.Ndim; Nc = size(u_sed,2);
            u = zeros(n,Nc,Ns); y = zeros(size(y_sed,1),Nc,Ns);
            Ncal=0;
            for k=1:Ns
                if k==1 && obj.SedSav
                    u(:,:,k)=u_sed; y(:,:,k)=y_sed; continue
                end
                nu=L*randn(n,Nc);
                threshold=y_sed(1,:)+log(rand(1,Nc));
                angle=2*pi*rand(1,Nc);
                lo=angle-2*pi; hi=angle;
                active=1:Nc;
                for attempt=1:obj.Ntry
                    proposed=mu+(u_sed(:,active)-mu).*cos(angle(active))+nu(:,active).*sin(angle(active));
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
