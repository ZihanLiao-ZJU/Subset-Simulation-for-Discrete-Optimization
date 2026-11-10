classdef MIS_Rel
    % Reliability MIS in the ISD's integration space. Default target is
    % I(g(u)<0)*p(u); projected views may supply EvlTarget(j) in their measure.
    properties
        TF
        Nite
        Pf
        logPf
        intBay
        X
        Y       % normalized: row 1 is logLi-logC_i
        Yraw    % raw INTEGRATION records: [logLi; logpi; domain payload]
        SamplingX
        SamplingYraw
        Context % same integration and target contract as MIS_Bay
        G
        logC
        W_res   % unscaled log ratios for sampled proposal transitions
        Pi_log
        Wmis
        Alpha
        Initial % external initial normalizer and its common-factor uncertainty
        ChainLengths % optional chain-major ragged layouts; empty means rectangular
        IndependentSamples = [] % conditional iid levels; [] retains legacy initial-only rule
        NormalizerInfo = struct('known',false(0,1),'exactLogC',[],'jacobian',[]) % normalizer dependency
    end
    methods
        function obj = MIS_Rel(in)
            obj.Context=MISContext(in,'reliability');
            obj.Nite=obj.Context.Nite;
            obj.ChainLengths=obj.Context.ChainLengths;
            obj.IndependentSamples=obj.Context.IndependentSamples;
            obj.Initial=obj.Context.Initial;
            obj.TF = obj.Context.TF;
            obj.SamplingX=obj.Context.SamplingX; obj.SamplingYraw=obj.Context.SamplingYraw;
            obj.intBay=obj.Context.IntBay; obj.X=obj.Context.X; obj.Yraw=obj.Context.Yraw;
            obj.G=obj.Context.G;
            [obj.logC,obj.Y,obj.W_res,obj.NormalizerInfo] = obj.normalizeY(obj.X,obj.Yraw,obj.G);
            assert(~obj.Initial.estimated || ~any(obj.NormalizerInfo.known), ...
                'SeMIS:InitialDependence','Exact absolute normalizers cannot use shared estimated-initial-scale propagation.');
            [obj.Wmis,obj.Alpha] = WgtMis_Blc(obj.Context.Views,obj.X,obj.Y,obj.G,obj.logC);
            n = obj.Context.TargetCount;
            obj.Pf = zeros(n,1);
            obj.logPf = -inf(n,1);
            obj.Pi_log = cell(n,1);
            for j = 1:n
                w = obj.targetWeights(j);
                obj.Pi_log{j} = MisInt(obj.Wmis,w);
                obj.logPf(j) = logsum(obj.Pi_log{j},'all');
                obj.Pf(j) = exp(obj.logPf(j));
            end
        end

        function [CoV_lw,CoV_up] = estimateCov(obj)
            % Same three-term v1 construction as Bayesian evidence; only f differs.
            CoV_lw = zeros(obj.TF.Nfun,1);
            CoV_up = zeros(obj.TF.Nfun,1);
            for j = 1:obj.TF.Nfun
                w = obj.targetWeights(j);
                for k = 1:obj.Nite, w{k} = w{k}+obj.Wmis{k}; end
                [CoV_lw(j),CoV_up(j)] = MisCov(obj.logPf(j),w,obj.Alpha,obj.W_res,obj.ChainLengths,obj.NormalizerInfo.jacobian,obj.IndependentSamples,obj.Context.SampleDesign);
                [CoV_lw(j),CoV_up(j)] = InitialNormalizerCov(CoV_lw(j),CoV_up(j),obj.Initial);
            end
        end

        function report = estimateUncertainty(obj,options)
            % Optional joint first-order estimate; legacy stopping is unchanged.
            if nargin<2, options=struct(); end
            reports=cell(obj.Context.TargetCount,1);
            [lo,hi]=obj.estimateCov;
            for j=1:obj.Context.TargetCount
                w=obj.targetWeights(j);
                for k=1:obj.Nite, w{k}=w{k}+obj.Wmis{k}; end
                r=MISUncertainty(obj.Context,obj.logPf(j),w,obj.Alpha, ...
                    obj.W_res,obj.logC,obj.NormalizerInfo.jacobian,options);
                r.legacyCovLower=lo(j); r.legacyCovUpper=hi(j);
                reports{j}=r;
            end
            if numel(reports)==1, report=reports{1}; else, report=reports; end
        end

        function [c_log,y_nor,w_res,info] = normalizeY(obj,x,y,g)
            [c_log,y_nor,w_res,info] = NormalizeMIS(obj.Context.Views,x,y,g,obj.Initial.logC);
        end
    end

    methods
        function w = targetWeights(obj,j)
            % Unbalanced log(target/q_i) values; shared by MIS, CoV and
            % diagnostic last-level IS, in the SAME integration measure.
            w=MISTargetWeights(obj.Context,obj.Y,j);
        end
    end
end
