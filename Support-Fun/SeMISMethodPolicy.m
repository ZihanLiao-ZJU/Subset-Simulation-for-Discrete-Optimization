classdef SeMISMethodPolicy
    % Optional, model-owned lifecycle hooks for the common SeMCS loop.
    % No hook (or an empty hook result) delegates to the existing engine path.
    %
    % SelectSeeds(context) returns parentIndices, candidateIndices, selection,
    % and plan. Indices refer to context.level.integration, flattened in
    % MATLAB column order. The policy obtains coordinates/cache from those
    % actual parents and performs the NEXT proposal's normal space lift.
    % A hook cannot substitute untracked seed coordinates or density rows.
    %
    % ShouldStop(context,phase) returns stop, converged and reason. A cap is
    % a stop, not convergence. Hooks must use cached context/model state;
    % stopping does not have a model-evaluation cost channel.
    methods (Static)
        function [result,model] = SelectSeeds(model,context)
            result=[];
            if ~ismethod(model,'SelectSeeds'), return; end
            assert(isstruct(context) && isscalar(context) && ...
                all(isfield(context,{'level','nextBudget','sampler'})), ...
                'SeMIS:MethodContext','Seed selection requires level, nextBudget and sampler.');
            validateattributes(context.nextBudget,{'numeric'},{'scalar','integer','positive','finite'});
            if isfield(context,'nextG')
                assert(isequaln(model.G,context.nextG),'SeMIS:MethodContext', ...
                    'The selection model must use the requested next ISD.');
            end
            level=context.level;
            assert(isstruct(level) && isscalar(level) && isfield(level,'integration') && ...
                isstruct(level.integration) && all(isfield(level.integration,{'x','y'})), ...
                'SeMIS:MethodContext','A frozen integration population is required.');
            request=model.SelectSeeds(context);
            if isempty(request), return; end
            assert(isstruct(request) && isscalar(request) && ...
                all(isfield(request,{'parentIndices','candidateIndices','selection','plan'})), ...
                'SeMIS:MethodSeeds','SelectSeeds must return actual parent indices, candidates, selection and plan.');
            assert(all(ismember(fieldnames(request), ...
                {'parentIndices','candidateIndices','selection','plan','metadata'})), ...
                'SeMIS:MethodSeeds','Seed coordinates and cached density rows are derived from parents, not supplied by a hook.');
            x=level.integration.x; y=level.integration.y;
            assert(isnumeric(x) && isreal(x) && isnumeric(y) && isreal(y) && ...
                size(y,1)>=3 && isequal(size(x,2:3),size(y,2:3)), ...
                'SeMIS:MethodContext','Integration coordinates and records must share a particle layout.');
            x=reshape(x,size(x,1),[]); y=reshape(y,size(y,1),[]);
            assert(all(isfinite(x),'all'),'SeMIS:MethodContext','Seed parents need finite integration coordinates.');
            n=size(x,2); idx=request.parentIndices; pool=request.candidateIndices;
            SeMISMethodPolicy.CheckIndices(idx,n,false);
            SeMISMethodPolicy.CheckIndices(pool,n,true);
            assert(all(ismember(idx,pool)),'SeMIS:MethodSeeds','Selected parents must belong to the declared candidate population.');
            selection=request.selection;
            assert(ischar(selection) && isrow(selection) && ...
                ismember(selection,{'weighted','survivors','custom'}), ...
                'SeMIS:MethodSeeds','Declare weighted, survivors or custom seed selection.');
            if strcmp(selection,'survivors')
                assert(numel(unique(idx))==numel(idx),'SeMIS:MethodSeeds', ...
                    'Survivor selection requires distinct parents.');
            end
            if isempty(idx)
                assert(isempty(request.plan),'SeMIS:MethodSeeds','An empty selection must not claim a chain plan.');
                result=struct('x',zeros(model.Ndim,0),'y',zeros(model.Nfun,0), ...
                    'cost',0,'plan',[],'seedInfo',[]);
                return
            end
            assert(~isempty(pool),'SeMIS:MethodSeeds','Nonempty seeds require a declared candidate population.');
            % Check NEXT-ISD support using cached payload only. These density
            % evaluators must not call EvlY/TF; no target cost is introduced.
            view=ISDIntegrationView(model);
            [logw,~,raw]=ResRat(view,x,y,model.G);
            logw=logw(:); raw=raw(:);
            assert(all(isfinite(y(1:2,idx)),'all') && all(isfinite(raw(idx))), ...
                'SeMIS:MethodSeeds','Selected parents must have positive density under the next ISD.');
            SeMISMethodPolicy.CheckDesign(model,x,y,logw,pool,selection);
            plan=SeMISMethodPolicy.ValidatePlan(request.plan,context.nextBudget);
            assert(plan.chains==numel(idx),'SeMIS:MethodSeeds','Chain count must equal the actual selected parent count.');
            if isfield(plan,'selection')
                assert(strcmp(plan.selection,selection),'SeMIS:MethodSeeds','The plan and seed design must agree.');
            end
            plan.selection=selection;
            sampler=context.sampler;
            retains=isprop(sampler,'SedSav') && sampler.SedSav;
            if isprop(sampler,'SeedMode'), retains=retains || ~strcmp(sampler.SeedMode,'move'); end
            assert(~strcmp(plan.keep,'last') || ~retains,'SeMIS:TransitionCount', ...
                'Terminal retention requires a sampler that does not count its input seed as a move.');
            seedInfo=struct('parentIndices',idx,'parentN',n,'selection',selection, ...
                'requestedSelection',selection,'candidateCount',numel(pool),'candidateIndices',pool, ...
                'integrationX',x(:,idx),'integrationY',y(:,idx), ...
                'integrationPreserved',false,'keep',plan.keep,'steps',plan.steps,'SeedMode','unknown');
            if isfield(request,'metadata')
                assert(isstruct(request.metadata) && isscalar(request.metadata), ...
                    'SeMIS:MethodSeeds','Custom seed metadata must be a scalar struct.');
                seedInfo.methodMetadata=request.metadata;
            end
            if isprop(sampler,'SeedMode'), seedInfo.SeedMode=char(sampler.SeedMode);
            elseif isprop(sampler,'SedSav')
                if sampler.SedSav, seedInfo.SeedMode='retain'; else, seedInfo.SeedMode='move'; end
            end
            space=ISDSpace(model); beforeG=model.G; dims=[model.Ndim,model.Nfun];
            seedInfo.integrationPreserved=strcmp(space.Kind,'identity') && ...
                isempty(space.Lift) && isempty(space.Project);
            if ismethod(model,'SeedTransportInvariant')
                invariant=model.SeedTransportInvariant();
                assert(islogical(invariant) && isscalar(invariant),'SeMIS:SeedTransportInvariant', ...
                    'SeedTransportInvariant must return a logical scalar.');
                seedInfo.integrationPreserved=seedInfo.integrationPreserved || invariant;
            end
            [xs,ys,model,cost]=space.LiftSamples(model,x(:,idx),y(:,idx));
            validateattributes(cost,{'numeric'},{'scalar','integer','nonnegative','finite'});
            assert(isequaln(beforeG,model.G) && isequal(dims,[model.Ndim,model.Nfun]), ...
                'SeMIS:LiftSchedule','Seed lifting must preserve the selected ISD and its sampling dimensions.');
            assert(isnumeric(xs) && isreal(xs) && isequal(size(xs),[model.Ndim,numel(idx)]) && ...
                all(isfinite(xs),'all') && isnumeric(ys) && isreal(ys) && ...
                isequal(size(ys),[model.Nfun,numel(idx)]) && all(isfinite(ys(1:2,:)),'all'), ...
                'SeMIS:MethodSeeds','Lifted seeds must have the current ISD layout and positive proposal density.');
            result=struct('x',xs,'y',ys,'cost',cost,'plan',plan,'seedInfo',seedInfo);
        end

        function decision = ShouldStop(model,context,phase)
            decision=[];
            if ~ismethod(model,'ShouldStop'), return; end
            assert(isstruct(context) && isscalar(context),'SeMIS:MethodContext','Stopping requires a scalar context.');
            assert(ischar(phase) && ismember(phase,{'before_seeds','after_seeds'}), ...
                'SeMIS:MethodPhase','Stopping phase must be before_seeds or after_seeds.');
            decision=model.ShouldStop(context,phase);
            if isempty(decision), return; end
            assert(isstruct(decision) && isscalar(decision) && ...
                all(isfield(decision,{'stop','converged','reason'})), ...
                'SeMIS:MethodStop','ShouldStop must return stop, converged and a canonical reason.');
            assert(islogical(decision.stop) && isscalar(decision.stop) && ...
                islogical(decision.converged) && isscalar(decision.converged) && ...
                ischar(decision.reason) && (isrow(decision.reason) || isempty(decision.reason)), ...
                'SeMIS:MethodStop','Stop/converged must be logical scalars and reason must be text.');
            if ~decision.stop
                valid=~decision.converged && isempty(decision.reason);
            elseif decision.converged
                valid=ismember(decision.reason,{'isd_complete','method_complete'});
            else
                valid=ismember(decision.reason,{'method_cap','method_stop'});
            end
            assert(valid,'SeMIS:MethodStop','The termination reason must distinguish convergence from a cap or other stop.');
        end

        function plan = ValidatePlan(plan,budget)
            % Same retention rules as SeMCS; custom only describes selection.
            assert(isstruct(plan) && isscalar(plan) && all(isfield(plan,{'chains','steps','keep'})), ...
                'SeMIS:SamplingPlan','A chain plan needs chains, steps and keep.');
            validateattributes(plan.chains,{'numeric'},{'scalar','integer','positive','finite'});
            validateattributes(plan.steps,{'numeric'},{'vector','integer','positive','finite'});
            assert(isscalar(plan.steps) || numel(plan.steps)==plan.chains, ...
                'SeMIS:SamplingPlan','A step vector needs one length per chain.');
            plan.keep=validatestring(plan.keep,{'all','last'});
            if ~isfield(plan,'budgetMode'), plan.budgetMode='exact'; end
            plan.budgetMode=validatestring(plan.budgetMode,{'exact','nominal'});
            if isfield(plan,'selection')
                plan.selection=validatestring(plan.selection,{'weighted','survivors','custom'});
            end
            retained=plan.chains;
            if strcmp(plan.keep,'all')
                if isscalar(plan.steps), retained=retained*plan.steps; else, retained=sum(plan.steps); end
            else
                assert(isscalar(plan.steps),'SeMIS:SamplingPlan','Terminal retention requires equal chain lengths.');
            end
            assert(retained==budget || strcmp(plan.budgetMode,'nominal'), ...
                'SeMIS:SamplingPlan','The chain plan must retain the requested exact budget.');
            plan.retained=retained;
        end
    end
    methods (Static, Access=private)
        function CheckDesign(model,x,y,logw,pool,selection)
            % Standard labels are declarations of the documented sampling
            % law, not inferred from one realised draw. Check all observable
            % support/design facts; arbitrary laws must use custom instead.
            if strcmp(selection,'custom'), return; end
            if strcmp(selection,'weighted')
                candidates=find(exp(logw)>0).';
            else
                active=isfinite(logw);
                scale=y(1:2,:); scale=scale(isfinite(scale));
                tolerance=32*eps(max([1;abs(scale(:))]));
                assert(all(abs(logw(active))<=tolerance),'SeMIS:MethodSeeds', ...
                    'Survivor selection requires constant finite transition weights.');
                if ismethod(model,'SeedEligibility')
                    eligible=model.SeedEligibility(x,y);
                    assert(islogical(eligible) && isvector(eligible) && ...
                        numel(eligible)==numel(active) && ~any(eligible(:) & ~active), ...
                        'SeMIS:SeedEligibility','SeedEligibility must select current integration support.');
                    active=active & eligible(:);
                end
                candidates=find(active).';
            end
            assert(isequal(sort(pool),sort(candidates)),'SeMIS:MethodSeeds', ...
                'A standard seed design must declare its complete eligible population; otherwise use custom.');
        end
        function CheckIndices(idx,n,distinct)
            assert(isnumeric(idx) && isreal(idx) && (isrow(idx) || isempty(idx)) && ...
                all(isfinite(idx) & idx==fix(idx) & idx>=1 & idx<=n), ...
                'SeMIS:MethodSeeds','Parent/candidate indices must index the actual flattened integration population.');
            if distinct
                assert(numel(unique(idx))==numel(idx),'SeMIS:MethodSeeds','Candidate indices must be distinct.');
            end
        end
    end
end
