classdef SuS_Multi_opt
    % Multiobjective discrete subset simulation in Gaussian sampling space.
    % Frozen-reference nondominated rank defines a pointwise indicator ISD.
    % Objectives are minimized; constraints, when supplied, satisfy C <= 0.
    properties
        OBJ
        Dim
        Nobj = 0
        Ncst = 0
        Nfun = 3
        Nsam = 500
        Ncal = 0
        X
        Y
        G
        Xmin
        Xmax
        Unique = false
        PermutationInvariant = false
        Mapping = 'auto'
        CacheZ
        CacheF
        CacheC
    end
    properties (Dependent)
        Ndim
        FlgCvg
        TF
    end
    methods
        function obj = SuS_Multi_opt(OBJ, xmin, xmax, varargin)
            options = inputParser;
            options.addParameter('Unique', []);
            options.addParameter('PermutationInvariant', []);
            options.addParameter('Mapping', 'auto');
            options.addParameter('ObjectiveCount', []);
            options.addParameter('ConstraintCount', []);
            options.addParameter('Patience', 5);
            options.addParameter('MinimumLevels', 3);
            options.addParameter('ConstraintTolerance', 0);
            options.parse(varargin{:});
            o = options.Results;
            obj.OBJ = OBJ;
            obj.Xmin = xmin(:); obj.Xmax = xmax(:);
            obj.Dim = numel(obj.Xmin);
            obj.Unique = obj.option(o.Unique, obj.field(OBJ, 'Unique', false));
            obj.PermutationInvariant = obj.option(o.PermutationInvariant, ...
                obj.field(OBJ, 'PermutationInvariant', false));
            if ~((islogical(obj.PermutationInvariant) || isnumeric(obj.PermutationInvariant)) && ...
                    isscalar(obj.PermutationInvariant) && isreal(obj.PermutationInvariant) && ...
                    isfinite(obj.PermutationInvariant) && any(obj.PermutationInvariant == [0, 1]))
                error('SuS_MultiOptim:PermutationInvariant', 'PermutationInvariant must be a scalar logical.');
            end
            obj.Mapping = char(o.Mapping);
            if strcmpi(obj.Mapping, 'auto')
                obj.Mapping = obj.field(OBJ, 'Mapping', 'auto');
            end
            DiscreteMap(zeros(obj.Dim, 0), obj.Xmin, obj.Xmax, obj.Unique, obj.Mapping);
            if obj.PermutationInvariant && ...
                    (~all(obj.Xmin == obj.Xmin(1)) || ~all(obj.Xmax == obj.Xmax(1)))
                error('SuS_MultiOptim:PermutationBounds', ...
                    'Permutation-invariant caching requires common bounds.');
            end
            dim = obj.field(OBJ, 'Ndim', obj.Dim);
            if ~isequal(dim, obj.Dim)
                error('SuS_MultiOptim:Dimension', 'Objective dimension and integer bounds disagree.');
            end
            obj.Nobj = obj.option(o.ObjectiveCount, obj.field(OBJ, 'Nfun', 0));
            obj.Ncst = obj.option(o.ConstraintCount, obj.field(OBJ, 'Ncst', 0));
            validateattributes(obj.Nobj, {'numeric'}, {'scalar','integer','nonnegative','finite'});
            validateattributes(obj.Ncst, {'numeric'}, {'scalar','integer','nonnegative','finite'});
            if obj.Nobj == 1
                error('SuS_MultiOptim:Objectives', 'Use SuS_Discrete_Optim for one objective.');
            end
            validateattributes(o.Patience, {'numeric'}, {'scalar','integer','positive','finite'});
            validateattributes(o.MinimumLevels, {'numeric'}, {'scalar','integer','positive','finite'});
            validateattributes(o.ConstraintTolerance, {'numeric'}, {'scalar','nonnegative','finite'});
            obj.Nfun = 3 + obj.Nobj + obj.Ncst;
            obj.CacheZ = zeros(obj.Dim, 0);
            obj.CacheF = zeros(obj.Nobj, 0); obj.CacheC = zeros(obj.Ncst, 0);
            obj.G = struct('p', .1, 'L', eye(obj.Dim), ...
                'referenceF', zeros(obj.Nobj, 0), 'referenceRanks', zeros(1, 0), ...
                'referenceViolation', zeros(1, 0), 'rankThreshold', inf, ...
                'iteration', 0, 'stableLevels', 0, 'patience', o.Patience, ...
                'minimumLevels', o.MinimumLevels, 'constraintTolerance', o.ConstraintTolerance, ...
                'archiveZ', zeros(obj.Dim, 0), 'archiveF', zeros(obj.Nobj, 0), ...
                'archiveC', zeros(obj.Ncst, 0));
        end

        function [u, y, Ncal, obj] = SamGen(obj, n)
            obj.checkReference();
            u = randn(obj.Dim, n);
            obj = obj.EvlY(u);
            y = [obj.EvlLKF;obj.EvlPDF;obj.Y];
            if any(~isfinite(y(1,:))) || any(abs(y(1,:))>1e-10)
                error('SeMIS:InitialISD','The first optimization ISD must have zero log likelihood.');
            end
            Ncal = obj.Ncal;
        end

        function Z = u2ints(obj, u)
            Z = DiscreteMap(u, obj.Xmin, obj.Xmax, obj.Unique, obj.Mapping);
        end

        function obj = EvlY(obj, u)
            obj.checkReference();
            u = reshape(u, obj.Dim, []);
            validateattributes(u, {'numeric'}, {'real','finite'});
            Z = obj.u2ints(u);
            if obj.PermutationInvariant, Z = sort(Z, 1); end
            [keys, ia, group] = unique(Z.', 'rows', 'stable');
            Zrep = Z(:, ia);
            [hit, location] = ismember(keys, obj.CacheZ.', 'rows');
            misses = find(~hit);
            obj.Ncal = numel(misses);
            if ~isempty(misses)
                [Fnew, Cnew] = obj.evaluate(Zrep(:, misses));
                if obj.Nobj == 0
                    obj.Nobj = size(Fnew, 1);
                    obj.CacheF = zeros(obj.Nobj, 0);
                    % Empty reference fields need no resize: the initial
                    % objective metadata must not change the requested ISD G.
                end
                if size(Fnew, 1) ~= obj.Nobj
                    error('SuS_MultiOptim:ObjectiveShape', 'Objective count changed or disagrees with Nfun.');
                end
                if size(Cnew, 1) ~= obj.Ncst
                    error('SuS_MultiOptim:ConstraintShape', 'Constraints must agree with declared Ncst.');
                end
                first = size(obj.CacheZ, 2);
                obj.CacheZ = [obj.CacheZ, Zrep(:, misses)];
                obj.CacheF = [obj.CacheF, Fnew];
                obj.CacheC = [obj.CacheC, Cnew];
                location(misses) = first + (1:numel(misses));
            end
            obj.Nfun = 3 + obj.Nobj + obj.Ncst;
            F = obj.CacheF(:, location);
            C = obj.CacheC(:, location);
            obj.X = u;
            obj.Y = [-.5 * (obj.Dim * log(2*pi) + sum(u.^2, 1)); F(:, group); C(:, group)];
        end

        function value = EvlLKF(obj)
            F = obj.Y(2:1+obj.Nobj, :);
            C = obj.Y(2+obj.Nobj:end, :);
            value = log(double(obj.accept(F, C)));
        end

        function value = EvlPDF(obj)
            value = obj.Y(1, :);
        end

        function [y, g, obj, Ncal] = UpdObj(obj, ~, y, p)
            if nargin < 4, p = obj.G.p; end
            validateattributes(p, {'numeric'}, {'scalar','>',0,'<',1});
            yflat = reshape(y, obj.Nfun, []);
            F = yflat(4:3+obj.Nobj, :);
            C = yflat(4+obj.Nobj:end, :);
            violation = sum(max(C - obj.G.constraintTolerance, 0), 1);
            ranks = NonDomRnk(F, violation);
            k = max(1, min(numel(ranks), round(numel(ranks) * p)));
            ordered = sort(ranks);
            retain = ranks <= ordered(k);
            % Whole rank fronts remain in support, including the selected seeds.
            [~, distinct] = unique([F(:, retain); violation(retain)].', 'rows', 'stable');
            pool = find(retain); pool = pool(distinct);
            g = obj.G;
            g.p = p; g.referenceF = F(:, pool);
            g.referenceRanks = ranks(pool); g.referenceViolation = violation(pool);
            g.rankThreshold = max(g.referenceRanks);
            g.iteration = g.iteration + 1;
            [Zarchive, Farchive, Carchive] = obj.paretoArchive();
            if isequal(sortrows(Farchive.'), sortrows(g.archiveF.'))
                g.stableLevels = g.stableLevels + 1;
            else
                g.stableLevels = 0;
            end
            g.archiveZ = Zarchive; g.archiveF = Farchive; g.archiveC = Carchive;
            obj.G = g;
            Ncal = 0;
        end

        function [xs, ys, chosen, candidates] = SedSlt(obj, x, y, n)
            % Rank/crowding seed choice is method-owned; kernels remain shared.
            x = reshape(x, obj.Dim, []); y = reshape(y, obj.Nfun, []);
            F = y(4:3+obj.Nobj, :); C = y(4+obj.Nobj:end, :);
            pool = find(obj.accept(F, C));
            candidates = pool;
            if isempty(pool)
                error('SuS_MultiOptim:Seeds', 'The selected ISD has no supported seed.');
            end
            nc = min(numel(pool), max(1, round(n * obj.G.p)));
            pool = pool(randperm(numel(pool)));
            violation = sum(max(C - obj.G.constraintTolerance, 0), 1);
            [~, uniqueValues] = unique([F(:, pool); violation(pool)].', 'rows', 'stable');
            distinct = pool(uniqueValues);
            if numel(distinct) >= nc
                chosen = obj.selectDiverse(F, violation, distinct, nc);
            else
                remainder = pool(~ismember(pool, distinct));
                chosen = [distinct, remainder(1:nc-numel(distinct))];
            end
            xs = x(:, chosen); ys = y(:, chosen);
            ys(1, :) = 0;
        end

        function request = SelectSeeds(obj, context)
            % The unified engine needs actual parents and the native chain plan.
            x = context.level.integration.x; y = context.level.integration.y;
            [~,~,parents,pool] = obj.SedSlt(x,y,context.nextBudget);
            nc = numel(parents);
            plan = struct('chains',nc,'steps',round(context.nextBudget/nc), ...
                'keep','all','budgetMode','nominal','selection','custom');
            request = struct('parentIndices',reshape(parents,1,[]), ...
                'candidateIndices',reshape(pool,1,[]),'selection','custom','plan',plan);
        end

        function [Z, F, C] = paretoArchive(obj)
            feasible = all(obj.CacheC <= obj.G.constraintTolerance, 1);
            Z = obj.CacheZ(:, feasible); F = obj.CacheF(:, feasible); C = obj.CacheC(:, feasible);
            if isempty(F), return; end
            retain = obj.firstFront(F);
            Z = Z(:, retain); F = F(:, retain); C = C(:, retain);
        end

        function d = get.Ndim(obj), d = obj.Dim; end
        function target = get.TF(obj)
            target = struct('Ndim',obj.Dim,'Domain','optimization');
        end
        function defaults = RunDefaults(~)
            defaults = struct('FlgCoV',false,'FlgCvgInt',true);
        end
        function yes = get.FlgCvg(obj)
            yes = obj.G.iteration >= obj.G.minimumLevels && ...
                ~isempty(obj.G.archiveF) && obj.G.stableLevels >= obj.G.patience;
        end
    end

    methods (Access=private)
        function checkReference(obj)
            factor = 1;
            if isfield(obj.G, 'L'), factor = obj.G.L;
            elseif isfield(obj.G, 'L_pi'), factor = obj.G.L_pi; end
            if ~isnumeric(factor) || ~isreal(factor) || ...
                    (~isequal(factor, 1) && ~isequal(factor, eye(obj.Dim)))
                error('SuS_MultiOptim:Reference', ...
                    'The discrete Gaussian pullback requires the standard-normal reference.');
            end
        end

        function mask = accept(obj, F, C)
            if isempty(obj.G.referenceF)
                mask = true(1, size(F, 2)); return;
            end
            violation = sum(max(C - obj.G.constraintTolerance, 0), 1);
            score = zeros(1, size(F, 2));
            for i = 1:size(obj.G.referenceF, 2)
                if obj.G.referenceViolation(i) == 0
                    better = all(obj.G.referenceF(:, i) <= F, 1) & ...
                        any(obj.G.referenceF(:, i) < F, 1);
                    dominated = violation > 0 | (violation == 0 & better);
                else
                    dominated = violation > obj.G.referenceViolation(i);
                end
                score(dominated) = max(score(dominated), obj.G.referenceRanks(i) + 1);
            end
            mask = score <= obj.G.rankThreshold;
        end

        function [F, C] = evaluate(obj, Z)
            provider = obj.OBJ;
            if isa(provider, 'function_handle')
                F = provider(Z);
            elseif isstruct(provider) && isfield(provider, 'EvlF')
                F = provider.EvlF(Z);
            elseif ismethod(provider, 'EvlF')
                F = provider.EvlF(Z);
            elseif isstruct(provider) && isfield(provider, 'EvlOBF')
                F = provider.EvlOBF(Z);
            elseif ismethod(provider, 'EvlOBF')
                F = provider.EvlOBF(Z);
            else
                error('SuS_MultiOptim:Objective', 'Provide a function handle, EvlF or EvlOBF.');
            end
            validateattributes(F, {'numeric'}, {'2d','real','finite'});
            if size(F, 1) < 2 || size(F, 2) ~= size(Z, 2)
                error('SuS_MultiOptim:ObjectiveShape', 'Return objectives-by-points with at least two objectives.');
            end
            F = double(F);
            C = zeros(0, size(Z, 2));
            if obj.Ncst > 0
                if isstruct(provider) && isfield(provider, 'EvlCST')
                    C = provider.EvlCST(Z);
                elseif ismethod(provider, 'EvlCST')
                    C = provider.EvlCST(Z);
                else
                    error('SuS_MultiOptim:Constraints', 'A positive Ncst requires EvlCST.');
                end
                validateattributes(C, {'numeric'}, {'2d','real','finite'});
                if ~isequal(size(C), [obj.Ncst, size(Z, 2)])
                    error('SuS_MultiOptim:ConstraintShape', 'Return Ncst-by-points constraints.');
                end
                C = double(C);
            end
        end
    end

    methods (Static, Access=private)
        function keep = firstFront(F)
            % Exact Pareto archive without a full cache-by-cache rank matrix.
            % Deduplicate first, then restore original cache order at the end.
            [~, first] = unique(F.', 'rows', 'stable');
            [ordered, order] = sortrows(F(:, first).');
            accepted = zeros(1, size(ordered, 1));
            count = 0;
            if size(F, 1) == 2
                best = inf;
                for k = 1:size(ordered, 1)
                    if ordered(k, 2) < best
                        count = count + 1; accepted(count) = k;
                        best = ordered(k, 2);
                    end
                end
            else
                for k = 1:size(ordered, 1)
                    if ~any(all(ordered(accepted(1:count), :) <= ordered(k, :), 2))
                        count = count + 1; accepted(count) = k;
                    end
                end
            end
            keep = sort(first(order(accepted(1:count))));
        end

        function value = field(provider, name, fallback)
            value = fallback;
            if isstruct(provider) && isfield(provider, name), value = provider.(name);
            elseif isobject(provider) && isprop(provider, name), value = provider.(name); end
        end
        function value = option(value, fallback)
            if isempty(value), value = fallback; end
        end
        function chosen = selectDiverse(F, violation, pool, n)
            rank = NonDomRnk(F(:, pool), violation(pool));
            chosen = [];
            for front = unique(rank)
                ids = pool(rank == front);
                remaining = n - numel(chosen);
                if numel(ids) <= remaining
                    chosen = [chosen, ids]; %#ok<AGROW>
                else
                    distance = zeros(1, numel(ids));
                    for objective = 1:size(F, 1)
                        [values, order] = sort(F(objective, ids));
                        span = values(end) - values(1);
                        if span <= 0, continue; end
                        distance(order([1, end])) = inf;
                        if numel(ids) > 2
                            distance(order(2:end-1)) = distance(order(2:end-1)) + ...
                                (values(3:end)-values(1:end-2)) / span;
                        end
                    end
                    [~, order] = sort(distance, 'descend');
                    chosen = [chosen, ids(order(1:remaining))]; %#ok<AGROW>
                end
                if numel(chosen) == n, break; end
            end
        end
    end
end
