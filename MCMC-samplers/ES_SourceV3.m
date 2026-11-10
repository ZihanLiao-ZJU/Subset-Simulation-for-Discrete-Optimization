classdef ES_SourceV3
    % Source Bayesian src v3 legacy coordinate-angle slice sampler.
    % Each coordinate draws its own angle on [0,pi]. This is NOT the
    % common-angle elliptical slice kernel implemented by ES_2pi.
    % Source mathematics, random-array layout and 20-try fallback are kept.
    % Only a standard Gaussian reference in U is supported. No aCS scale,
    % lambda carry, score adaptation or ragged chain lengths are applicable.
    properties (Constant)
        SedSav = false
        Ntry = 20
    end
    properties
        BayStc
    end
    methods
        function obj = ES_SourceV3(baystc)
            if nargin>0, obj.BayStc=baystc; end
        end

        function obj = ConfigureISDPolicy(obj,defaults)
            % SeMCS validates all legal policy names and values first. Both
            % kernels support the same move-only, non-aCS capability policy.
            checker=ES_2pi;
            checker.ConfigureISDPolicy(defaults);
        end

        function [u,y,Ncal,baystc] = SamGen(obj,varargin)
            % Seeded use moves each seed. The source unseeded overload
            % SamGen(Nc,Ns) includes its generated initial states in slot 1.
            if numel(varargin)~=2 && numel(varargin)~=3
                error('SeMIS:SamplerArguments','Use (Nc,Ns) or (u_sed,y_sed,Ns).');
            end
            Ns=varargin{end};
            if ~isscalar(Ns)
                error('SeMIS:RaggedSampler','ES_SourceV3 requires one common chain length.');
            end
            validateattributes(Ns,{'numeric'},{'real','integer','positive','finite'});
            baystc=obj.BayStc;
            [mu,L]=GaussianPrior(baystc);
            if any(mu~=0) || ~isequal(L,eye(baystc.Ndim))
                error('SeMIS:SourceV3Prior', ...
                    'The source coordinate-angle kernel requires a standard Gaussian reference in U.');
            end
            sedsav=obj.SedSav || nargin<=3;
            nTry=obj.Ntry;
            nDim=baystc.Ndim;
            nFun=baystc.Nfun;
            if nargin<=3
                Nc=varargin{1};
                validateattributes(Nc,{'numeric'},{'scalar','real','integer','positive','finite'});
                u_sed=randn(nDim,Nc);
                baystc=baystc.EvlY(u_sed);
                y_sed=[baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                Ncal=baystc.Ncal;
            else
                u_sed=varargin{1};
                y_sed=varargin{2};
                Nc=size(u_sed,2);
                if ~isequal(size(u_sed),[nDim,Nc]) || ~isequal(size(y_sed),[nFun,Nc])
                    error('SeMIS:SamplerLayout','Seed coordinates and records must be two-dimensional arrays with matching columns.');
                end
                Ncal=0;
            end

            u=zeros(nDim,Nc,Ns);
            y=zeros(nFun,Nc,Ns);
            du=randn(nDim,Nc,Ns);
            p=rand(1,Nc,Ns);

            for k=1:Ns
                if k==1 && sedsav
                    u(:,:,1)=u_sed;
                    y(:,:,1)=y_sed;
                else
                    Rho=sqrt(u_sed.^2+du(:,:,k).^2);
                    alpha_appro=acos(u_sed./Rho);
                    ind_sam=1:Nc;
                    Nrst=Nc;
                    alphamin=zeros(nDim,Nc);
                    alphamax=pi*ones(nDim,Nc);
                    for i=1:nTry
                        alpha=rand(nDim,Nrst).*(alphamax-alphamin)+alphamin;
                        u(:,ind_sam,k)=cos(alpha).*Rho(:,ind_sam);
                        baystc=baystc.EvlY(u(:,ind_sam,k));
                        y(:,ind_sam,k)=[baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                        Ncal=Ncal+baystc.Ncal;
                        ind_rej=exp(y(1,ind_sam,k)-y_sed(1,ind_sam))<p(1,ind_sam,k);
                        Nrst=sum(ind_rej);
                        if Nrst>0
                            ind_sam=ind_sam(ind_rej);
                            u(:,ind_sam,k)=u_sed(:,ind_sam);
                            y(:,ind_sam,k)=y_sed(:,ind_sam);
                            alpha=alpha(:,ind_rej);
                            alphamin=alphamin(:,ind_rej);
                            alphamax=alphamax(:,ind_rej);
                            ind_alpha=alpha>alpha_appro(:,ind_sam);
                            alphamin(~ind_alpha)=alpha(~ind_alpha);
                            alphamax(ind_alpha)=alpha(ind_alpha);
                        else
                            break
                        end
                    end
                    u_sed=u(:,:,k);
                    y_sed=y(:,:,k);
                end
            end
            % Preserve the source return object (last evaluated proposal
            % batch). SeMCS owns assembly of retained X/Y after sampling.
        end
    end
end
