classdef SuS_Discrete_Optim
    %SUS_DISCRETE_OPTIM  Custom ISD (IntBay) for discrete optimisation with SeMCS.
    %
    %   SuS_Discrete_Optim is the custom intermediate sampling distribution (ISD,
    %   i.e. the framework IntBay) for discrete optimisation: the object stored
    %   in sampler.BayStc and returned as out.intBay.
    %
    %   SeMCS samples u in R^Dim from a standard-normal prior and accepts /
    %   rejects proposals on the level set {F(theta(u)) < threshold}.  The
    %   mapping from u to a discrete point is
    %
    %       U  = Phi(u)                       uniform in (0,1)^Dim
    %       Xc = xmin + (xmax-xmin).*U        continuous point in the box
    %       theta = u2ints(Xc)                discrete point (see 'Unique')
    %
    %   'Unique' (default true): all coordinates must be pairwise distinct
    %   integers (a subset / combination; cont2uniqint is used).
    %   'Unique' false: each coordinate is rounded independently inside its
    %   own [xmin(i), xmax(i)] bounds (a general integer program).
    %
    %   Objective interface (minimisation, F(theta) -> scalar values):
    %       function handle        F = OBJ(theta)
    %       object/struct with    F = OBJ.EvlF(theta)
    %       object with           F = OBJ.EvlOBF(theta)
    %   where theta is [Dim x K] integers and F is [1 x K].
    %
    %   SeMCS protocol implemented here:
    %       SamGen(N)        [u,y,Ncal,model] for the initial Gaussian level
    %       EvlY(u)          objective evaluation, Y = [log_phi; F]
    %       EvlLKF           log 1{F < G.parameters}
    %       EvlPDF           log phi_Dim(u)
    %       UpdObj(u,y,p)    new threshold = p-quantile of F
    %       FlgCvg           convergence flag (threshold collapse / no progress)
    %       Dim, Nfun, Ncal, X, Y, G, Nsam
    %
    %   Usage:
    %       model = SuS_Discrete_Optim(OBJ, xmin, xmax);
    %       model.Nsam = 150; model.G.p = 0.1;
    %       sampler = aCS; sampler.BayStc = model;
    %       simulation = SeMCS(sampler); simulation.MaxIte = 30;
    %       raw = simulation.RunIte;
    %
    %   See also SUS_MULTI_OPT, SEMCS, DISCRETEMAP.

    properties
        OBJ                 % objective provider (handle / object / struct)
        Nfun = 4            % packed rows: [log L; log pi; log_phi; F]
        Dim                 % number of integer variables
        Ncal = 0            % number of NEW objective evaluations in last EvlY
        X                   % last evaluated samples in u-space
        Y                   % [2 x N]: [log_phi; F]
        G                   % intermediate state struct (parameters,p,L,...)
        Nsam = 150          % samples per level (read by SeMCS constructor)
        ThresholdRule = 'hazen'     % 'hazen' (quantile) or 'midpoint'
        Unique = true               % pairwise-distinct integers (subset mode)
        PermutationInvariant = true % objective depends only on the multiset
        Mapping = 'auto'            % 'auto' | 'quantile' | 'rank' | 'uniform'
        Xmin                % per-dimension lower integer bound
        Xmax                % per-dimension upper integer bound
        CacheZ              % cached keys (sorted if permutation invariant)
        CacheF = double.empty(1,0)   % cached objective values [1 x K]
    end

    properties (Dependent)
        FlgCvg
        Ndim        % alias of Dim (framework protocol name)
        TF          % dimension metadata for the unified project engine
    end

    methods
        function obj = SuS_Discrete_Optim(OBJ, xmin, xmax, varargin)
            p = inputParser;
            p.addParameter('Unique', [], @(x)isempty(x)||((islogical(x)||isnumeric(x))&&isscalar(x)&&isreal(x)&&isfinite(x)&&(x==0||x==1)));
            p.addParameter('PermutationInvariant', [], @(x)isempty(x)||((islogical(x)||isnumeric(x))&&isscalar(x)&&isreal(x)&&isfinite(x)&&(x==0||x==1)));
            p.addParameter('ThresholdRule', 'hazen', @(x)ischar(x)||isstring(x));
            p.addParameter('Mapping', 'auto', @(x)ischar(x)||isstring(x));
            p.addParameter('StopRule','legacy',@(x)ischar(x)||isstring(x));
            p.addParameter('Patience',5,@(x)isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=1&&x==fix(x));
            p.addParameter('MinimumLevels',3,@(x)isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=1&&x==fix(x));
            p.parse(varargin{:});
            opt = p.Results;

            obj.OBJ = OBJ;
            validateattributes(xmin, {'numeric'}, {'real','vector','finite','integer','nonempty'});
            validateattributes(xmax, {'numeric'}, {'real','vector','finite','integer','nonempty'});
            obj.Xmin = double(xmin(:).');
            obj.Xmax = double(xmax(:).');
            if numel(obj.Xmin) ~= numel(obj.Xmax)
                error('SuS_DiscreteOptim:Bounds', 'xmin and xmax must have the same length.');
            end
            obj.Dim = numel(obj.Xmin);

            if isempty(opt.Unique)
                obj.Unique = obj.objField_(OBJ, 'Unique', true);
            else
                obj.Unique = logical(opt.Unique);
            end
            assert(isscalar(obj.Unique) && (islogical(obj.Unique)||isnumeric(obj.Unique)) && ...
                isreal(obj.Unique) && isfinite(obj.Unique) && any(obj.Unique==[0,1]), ...
                'SuS_DiscreteOptim:Unique', 'Unique must be a scalar logical.');
            commonBounds=all(obj.Xmin==obj.Xmin(1)) && all(obj.Xmax==obj.Xmax(1));
            permutationBounds = obj.Unique && all(obj.Xmin == 1) && all(obj.Xmax == obj.Dim);
            if isempty(opt.PermutationInvariant)
                obj.PermutationInvariant = obj.objField_(OBJ, 'PermutationInvariant', obj.Unique && commonBounds && ~permutationBounds);
            else
                obj.PermutationInvariant = logical(opt.PermutationInvariant);
            end
            obj.ThresholdRule = lower(char(opt.ThresholdRule));
            assert(any(strcmp(obj.ThresholdRule, {'hazen','midpoint'})), ...
                'SuS_DiscreteOptim:ThresholdRule', 'ThresholdRule must be hazen or midpoint.');
            if strcmpi(char(opt.Mapping), 'auto')
                m = obj.objField_(OBJ, 'Mapping', 'auto');
                obj.Mapping = lower(char(m));
            else
                obj.Mapping = lower(char(opt.Mapping));
            end
            assert(any(strcmp(obj.Mapping, {'auto','quantile','rank','uniform'})), ...
                'SuS_DiscreteOptim:Mapping', 'Mapping must be auto, quantile, rank or uniform.');
            assert(isscalar(obj.PermutationInvariant) && ...
                (islogical(obj.PermutationInvariant)||isnumeric(obj.PermutationInvariant)) && ...
                isreal(obj.PermutationInvariant) && isfinite(obj.PermutationInvariant) && ...
                any(obj.PermutationInvariant==[0,1]), ...
                'SuS_DiscreteOptim:PermutationInvariant', 'PermutationInvariant must be a scalar logical.');
            obj.Unique=logical(obj.Unique);
            obj.PermutationInvariant=logical(obj.PermutationInvariant);
            if obj.PermutationInvariant && ~commonBounds
                error('SuS_DiscreteOptim:PermutationInvariant', ...
                    'Canonical sorting requires common bounds in all coordinates.');
            end
            DiscreteMap(zeros(obj.Dim,0),obj.Xmin,obj.Xmax,obj.Unique,obj.Mapping);

            obj.CacheZ = zeros(obj.Dim, 0);
            obj.G = struct('parameters', inf, 'p', 0.1, 'L', eye(obj.Dim), ...
                'Ite', 0, 'sedF', double.empty(1,0), ...
                'NoProgress', false, 'prev', inf);
            stopRule=validatestring(opt.StopRule,{'legacy','stability'});
            if strcmp(stopRule,'stability')
                obj.G.StopRule=stopRule;
                obj.G.Patience=opt.Patience;
                obj.G.MinimumLevels=opt.MinimumLevels;
                obj.G.best=inf;
                obj.G.stableLevels=0;
            end
        end

        function [u, y, Ncal, obj] = SamGen(obj, N)
            % Prior sampling in standard-normal space (framework contract).
            factor=1;
            if isfield(obj.G,'L'),factor=obj.G.L;
            elseif isfield(obj.G,'L_pi'),factor=obj.G.L_pi;end
            if ~isequal(factor,1) && ~isequal(factor,eye(obj.Dim))
                error('SuS_DiscreteOptim:Reference', ...
                    'Discrete optimization requires the standard-normal reference.');
            end
            u = randn(obj.Dim, N);
            obj = obj.EvlY(u);
            y = [obj.EvlLKF;obj.EvlPDF;obj.Y];
            if any(~isfinite(y(1,:))) || any(abs(y(1,:))>1e-10)
                error('SeMIS:InitialISD','The first optimization ISD must have zero log likelihood.');
            end
            Ncal = obj.Ncal;
        end

        function obj = EvlY(obj, u)
            % Map u -> discrete point, score with cache, pack Y = [log_phi; F].
            factor=1;
            if isfield(obj.G,'L'),factor=obj.G.L;
            elseif isfield(obj.G,'L_pi'),factor=obj.G.L_pi;end
            if ~isequal(factor,1) && ~isequal(factor,eye(obj.Dim))
                error('SuS_DiscreteOptim:Reference', ...
                    'Discrete optimization requires the standard-normal reference.');
            end
            if ndims(u) > 2
                u = reshape(u, obj.Dim, []);
            end
            obj.X = u;
            n = size(u, 2);
            log_phi = -.5 * (obj.Dim * log(2*pi) + sum(u.^2, 1));

            Z = obj.u2ints(u);                      % [Dim x n] integers
            Zkey = obj.cacheKey_(Z);
            [~, ia, ic] = unique(Zkey.', 'rows', 'stable');
            Zrep = Z(:, ia);
            Nuniq = size(Zrep, 2);

            [in_cache, locb] = ismember(Zkey(:, ia).', obj.CacheZ.', 'rows');
            in_cache = in_cache.';
            locb = locb.';
            Fu_rep = nan(1, Nuniq);
            Fu_rep(in_cache) = obj.CacheF(locb(in_cache));

            miss = find(~in_cache);
            if isempty(miss)
                obj.Ncal = 0;
            else
                Zeval = Zrep(:, miss);
                if obj.PermutationInvariant
                    Zeval = sort(Zeval, 1);   % canonical representative
                end
                Fnew = obj.evalObjective_(Zeval);
                Fu_rep(miss) = Fnew;
                obj.CacheZ = [obj.CacheZ, Zkey(:, ia(miss))];
                obj.CacheF = [obj.CacheF, Fnew];
                obj.Ncal = numel(miss);
            end

            F = Fu_rep(ic);
            F = reshape(F, 1, n);
            obj.Y = [log_phi; F];
        end

        function log_L = EvlLKF(obj)
            % Log indicator likelihood: log 1{F < G.parameters}.
            F = obj.Y(2, :);
            F = F(:).';
            if obj.G.parameters == inf
                log_L = zeros(1, numel(F));
            else
                log_L = log(double(F < obj.G.parameters));
            end
        end

        function log_pi = EvlPDF(obj)
            % Standard-normal prior density in u-space.
            X = obj.X;
            if ndims(X) > 2
                X = reshape(X, obj.Dim, []);
            end
            log_pi = -.5 * (obj.Dim * log(2*pi) + sum(X.^2, 1));
        end

        function [y, g, obj, N_cal] = UpdObj(obj, u, y, p)
            if nargin < 4, p = obj.G.p; end
            % Update the level threshold from the current objective values.
            % The packed sampled y is returned unchanged: the generic
            % SedSlt/ResRat compares the new ISF against y(1:2,:).
            F = y(4, :);
            F = F(:).';
            u = reshape(u, obj.Dim, []);
            Yc = y(3:end, :);
            obj.X = u;
            obj.Y = reshape(Yc, size(Yc, 1), numel(F));

            g = obj.G;
            prev = g.parameters;
            switch obj.ThresholdRule
                case 'midpoint'
                    fs = sort(F, 'ascend');
                    k = max(1, min(numel(fs), round(numel(fs) * p)));
                    if k >= numel(fs)
                        l = fs(end);
                    else
                        l = (fs(k) + fs(k+1)) / 2;
                    end
                otherwise
                    l = quantile(F, p);
            end
            if isfield(g,'StopRule') && strcmp(g.StopRule,'stability')
                best=min(obj.CacheF);
                if best<g.best,g.stableLevels=0;else,g.stableLevels=g.stableLevels+1;end
                g.best=best;
                % A discrete quantile can coincide with every retained value.
                % Keep those ties inside the strict event instead of emptying it.
                if ~any(F<l),l=min(prev,l+eps(max(1,abs(l))));end
            end
            g.Ite = g.Ite + 1;
            g.p = p;
            g.parameters = l;
            g.sedF = F(F < l);
            g.NoProgress = isfinite(prev) && ...
                (l >= prev - 1e-12 * max(1, abs(prev)));
            g.prev = l;
            obj.G = g;

            g = obj.G; %#ok<NASGU>
            N_cal = 0;
        end

        function d = get.Ndim(obj)
            d = obj.Dim;
        end

        function target = get.TF(obj)
            target = struct('Ndim',obj.Dim,'Domain','optimization');
        end

        function defaults = RunDefaults(~)
            % Optimization has no Bayesian evidence/CoV stopping target.
            defaults = struct('FlgCoV',false,'FlgCvgInt',true);
        end

        function request = SelectSeeds(obj, context)
            % Actual binary-event parent selection, with nominal chain rounding.
            % Parent indices avoid dependence on another toolbox's SedSlt file.
            y = reshape(context.level.integration.y,obj.Nfun,[]);
            if obj.G.parameters == inf, pool = 1:size(y,2);
            else, pool = find(y(4,:) < obj.G.parameters); end
            n = context.nextBudget;
            if isempty(pool)
                parents = zeros(1,0); plan = [];
            else
                counts = round(n ./ round(n ./ (1:n)));
                nc = counts(find(counts<=numel(pool),1,'last'));
                parents = pool(randperm(numel(pool),nc));
                plan = struct('chains',nc,'steps',round(n/nc),'keep','all', ...
                    'budgetMode','nominal','selection','custom');
            end
            request = struct('parentIndices',reshape(parents,1,[]), ...
                'candidateIndices',reshape(pool,1,[]),'selection','custom','plan',plan);
        end

        function decision = ShouldStop(obj, context, phase)
            decision = [];
            % Explicitly disabling intrinsic convergence must remain effective.
            if ~context.runtimeG.run.flag_isd, return; end
            stop = strcmp(phase,'after_seeds') && obj.FlgCvg;
            reason = '';
            if stop, reason = 'method_complete'; end
            decision = struct('stop',logical(stop),'converged',logical(stop),'reason',reason);
        end

        function obj = set.Ndim(obj, d)
            obj.Dim = d;
        end

        function Flg = get.FlgCvg(obj)
            % Converged when the new level no longer shrinks (discrete
            % objective) or the surviving seed values collapse.
            if isfield(obj.G,'StopRule') && strcmp(obj.G.StopRule,'stability')
                Flg=obj.G.Ite>=obj.G.MinimumLevels && obj.G.stableLevels>=obj.G.Patience;
                return;
            end
            if isinf(obj.G.parameters)
                Flg = false;
                return;
            end
            sed = obj.G.sedF;
            if isempty(sed)
                Flg = true;
                return;
            end
            tol = 1e-12 * max(1, max(abs(sed)));
            Flg = (max(sed) - min(sed) <= tol) || obj.G.NoProgress;
        end

        function Z = u2ints(obj, u)
            % The mathematical map is shared with the multiobjective ISD.
            Z = DiscreteMap(u,obj.Xmin,obj.Xmax,obj.Unique,obj.Mapping);
        end

        function F = evalObjective_(obj, Z)
            % Dispatch to the user objective: handle / EvlF / legacy EvlOBF.
            C = obj.OBJ;
            if isa(C, 'function_handle')
                F = C(Z);
            elseif ismethod(C, 'EvlF')
                F = C.EvlF(Z);
            elseif ismethod(C, 'EvlOBF')          % OSP back-compatibility
                F = C.EvlOBF(Z);
            elseif isstruct(C) && isfield(C, 'EvlF')
                F = C.EvlF(Z);
            elseif isstruct(C) && isfield(C, 'EvlOBF')
                F = C.EvlOBF(Z);
            else
                error('SuS_DiscreteOptim:Objective', ...
                    'OBJ must be a function handle or provide EvlF (or legacy EvlOBF).');
            end
            if ~(isnumeric(F)||islogical(F)) || ~isreal(F) || ~isvector(F) || ...
                    numel(F)~=size(Z,2) || any(~isfinite(F(:)))
                error('SuS_DiscreteOptim:Objective', ...
                    'The objective must return one finite real value per design column.');
            end
            F = double(reshape(F, 1, []));
        end

        function K = cacheKey_(obj, Z)
            if obj.PermutationInvariant
                K = sort(Z, 1);
            else
                K = Z;
            end
        end

        function value = objField_(~, C, name, default)
            % Read a property/field from the user objective if present.
            value = default;
            if isobject(C) && isprop(C, name)
                value = C.(name);
            elseif isstruct(C) && isfield(C, name)
                value = C.(name);
            end
        end
    end

end
