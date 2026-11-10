classdef MIS_Bay
    % Bayesian MIS in the declared integration measure. Z is log evidence.
    % EvlTarget defines custom targets; legacy payload is [logL; logprior].
    properties
        X
        Y       % normalized proposal records: row 1 is logLi-logC_i
        Yraw    % raw INTEGRATION cache [logLi; logpi; declared target payload]
        SamplingX
        SamplingYraw
        Context % common integration history and frozen per-level views
        NumIte
        Wmis
        Wpos    % log posterior/proposal density ratios (not particle masses)
        Z       % log evidence
        logZ    % log normalizers of the NumIte actually sampled proposals
        Wres    % raw log ratios between successive sampled proposals
        G
        IntBay
        Alpha   % cross-level log balance weights for CoV derivatives
        Initial % external initial normalizer and its common-factor uncertainty
        ChainLengths % optional chain-major ragged layouts
        IndependentSamples = [] % conditional iid levels; [] retains legacy initial-only rule
        NormalizerInfo = struct('known',false(0,1),'exactLogC',[],'jacobian',[]) % normalizer dependency
    end
    methods
        function obj = MIS_Bay(in)
            obj.Context=MISContext(in,'bayesian');
            obj.X=obj.Context.X; obj.Yraw=obj.Context.Yraw;
            obj.SamplingX=obj.Context.SamplingX; obj.SamplingYraw=obj.Context.SamplingYraw;
            obj.ChainLengths=obj.Context.ChainLengths;
            obj.IndependentSamples=obj.Context.IndependentSamples;
            obj.Initial=obj.Context.Initial; obj.NumIte=obj.Context.Nite;
            obj.G=obj.Context.G; obj.IntBay=obj.Context.IntBay;
            [obj.logZ,obj.Y,obj.Wres,obj.NormalizerInfo] = NormalizeMIS(obj.Context.Views,obj.X,obj.Yraw,obj.G,obj.Initial.logC);
            assert(~obj.Initial.estimated || ~any(obj.NormalizerInfo.known), ...
                'SeMIS:InitialDependence','Exact absolute normalizers cannot use shared estimated-initial-scale propagation.');
            [obj.Wmis,obj.Alpha] = WgtMis_Blc(obj.Context.Views,obj.X,obj.Y,obj.G,obj.logZ);
            obj.Wpos = obj.targetWeights;
            obj.Z = logsum(MisInt(obj.Wmis,obj.Wpos),'all');
            if isfinite(obj.Z)
                for k = 1:obj.NumIte, obj.Wpos{k} = obj.Wpos{k}-obj.Z; end
            else
                % Zero evidence has undefined posterior; preserve zero mass so
                % PosSim/KLDEst raise a clear error rather than propagate NaNs.
                for k = 1:obj.NumIte, obj.Wpos{k}(:) = -inf; end
            end
        end

        function [CoV_lw,CoV_up] = estimateCov(obj)
            w = obj.targetWeights;
            for k = 1:obj.NumIte, w{k} = w{k}+obj.Wmis{k}; end
            [CoV_lw,CoV_up] = MisCov(obj.Z,w,obj.Alpha,obj.Wres,obj.ChainLengths,obj.NormalizerInfo.jacobian,obj.IndependentSamples,obj.Context.SampleDesign);
            [CoV_lw,CoV_up] = InitialNormalizerCov(CoV_lw,CoV_up,obj.Initial);
        end

        function report = estimateUncertainty(obj,options)
            % Optional joint first-order estimate; legacy stopping is unchanged.
            if nargin<2, options=struct(); end
            w=obj.targetWeights;
            for k=1:obj.NumIte, w{k}=w{k}+obj.Wmis{k}; end
            report=MISUncertainty(obj.Context,obj.Z,w,obj.Alpha,obj.Wres, ...
                obj.logZ,obj.NormalizerInfo.jacobian,options);
            [report.legacyCovLower,report.legacyCovUpper]=obj.estimateCov;
        end

        function w = targetWeights(obj)
            % Unnormalized evidence target, shared by integration and CoV.
            w=MISTargetWeights(obj.Context,obj.Y,1);
        end

        function [u_pos,y_pos,u_all,y_all] = PosSim(obj,Npos,space)
            % Posterior masses are computed in integration coordinates.
            % Explicit posterior recovery uses the TARGET conditional;
            % proposal Lift is never used. Outputs 3/4 stay integration data.
            if nargin<3, space='integration'; end
            space=validatestring(space,{'integration','posterior','output'});
            [u_all,y_all] = VabCmb(obj.X,obj.Yraw);
            logMass = MisParticleWeights(obj.Wmis,obj.Wpos);
            Nall = numel(logMass);
            if nargin<2 || isempty(Npos)
                Npos = min(Nall,max(1,round(exp(-logsum(2*logMass,'all')))));
            end
            validateattributes(Npos,{'numeric'},{'scalar','integer','positive','finite'});
            indices = randsample(Nall,Npos,true,exp(logMass));
            u_pos = u_all(:,indices);
            y_pos = y_all(:,indices);
            if ~strcmp(space,'integration')
                sizes=cellfun(@(y)prod(size(y,2:3)),obj.Yraw);
                origins=repelem(1:obj.NumIte,sizes(:).');
                origins=origins(indices);
                inputU=u_pos; inputY=y_pos;
                allocated=false;
                for k=1:obj.NumIte
                    selected=find(origins==k);
                    if isempty(selected), continue; end
                    adapter=obj.Context.Spaces{k}; model=obj.Context.Models{k};
                    if strcmp(space,'posterior')
                        [u,y]=adapter.PosteriorSamples(model,inputU(:,selected),inputY(:,selected));
                        dimension=obj.Context.Descriptors{k}.sampling.dimension;
                    else
                        [u,y]=adapter.OutputSamples(model,inputU(:,selected),inputY(:,selected));
                        dimension=obj.Context.Descriptors{k}.output.dimension;
                    end
                    assert(isnumeric(u) && isreal(u) && isequal(size(u),[dimension,numel(selected)]) && ...
                        all(isfinite(u),'all') && isnumeric(y) && isreal(y) && ismatrix(y) && ...
                        size(y,2)==numel(selected), ...
                        'SeMIS:PosteriorLayout','Posterior/output mappings must preserve particle columns and return valid coordinates and cache.');
                    if ~allocated
                        u_pos=zeros(dimension,Npos); y_pos=zeros(size(y,1),Npos);
                        allocated=true;
                    end
                    assert(size(y_pos,1)==size(y,1),'SeMIS:PosteriorLayout', ...
                        'Posterior/output payload dimensions must agree across sampled levels.');
                    u_pos(:,selected)=u; y_pos(:,selected)=y;
                end
            end
        end

        function [ccdf,dcdf] = CCDF(obj,L)
            % Legacy outputs are LOG probabilities: prior P(logL>=L_j), and
            % posterior interval masses [L_(j-1),L_j), starting at -Inf.
            obj.requireLikelihoodDiagnostics(true);
            assert(isreal(L) && ~any(isnan(L(:))) && all(diff(L(:))>=0), ...
                'SeMIS:Thresholds','Likelihood thresholds must be sorted real values.');
            assert(isfinite(obj.Z),'SeMIS:ZeroPosterior','Posterior requires positive evidence.');
            ccdf = zeros(size(L));
            dcdf = zeros(size(L));
            for j = 1:numel(L)
                ccdf(j) = logsum(MisInt(obj.Wmis,IndPriWgt(obj.Y,L(j))),'all');
                if j==1, lower = -inf; else, lower = L(j-1); end
                dcdf(j) = logsum(MisInt(obj.Wmis,IndPosWgt(obj.Y,obj.Z,lower,L(j))),'all');
            end
        end

        function [Ip,BMD] = KLDEst(obj)
            % E_post[logL-logZ] is signed pointwise. Accumulate real moments,
            % avoiding the old complex logarithm of a possibly negative term.
            obj.requireLikelihoodDiagnostics(false);
            logMass = MisParticleWeights(obj.Wmis,obj.Wpos);
            [~,y] = VabCmb(obj.X,obj.Yraw);
            keep = isfinite(logMass);
            p = exp(logMass(keep));
            d = y(3,keep)-obj.Z;
            Ip = sum(p.*d);
            BMD = 2*sum(p.*(d-Ip).^2);
        end

        function Nunq = UnqSam(obj)
            Nunq = zeros(obj.NumIte,1);
            for k = 1:obj.NumIte
                u = obj.X{k}(1:end,:);
                Nunq(k) = size(unique(u.','rows','stable'),1);
            end
        end
    end
    methods (Access=private)
        function requireLikelihoodDiagnostics(obj,priorEvents)
            for k=1:obj.NumIte
                descriptor=obj.Context.Descriptors{k};
                assert(descriptor.likelihoodDiagnostics && size(obj.Yraw{k},1)>=4 && ...
                    (~priorEvents || strcmp(descriptor.kind,'identity')), ...
                    'SeMIS:UnsupportedObservable', ...
                    'This integration space does not declare the cached likelihood/prior observable required by this diagnostic.');
            end
        end
    end
end
