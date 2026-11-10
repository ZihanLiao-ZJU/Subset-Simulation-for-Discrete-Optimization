classdef SeMCS
    % Sequential sampling engine for q_i = exp(logLi_i + logpi_i) / C_i.
    % Sampling and integration may use different spaces. Optional ISD hooks
    % ProjectSamples/IntegrationView/LiftSamples connect them; ordinary ISDs
    % use the identity map. The sampler always receives sampling-space data.
    properties (Dependent)
        % All public configuration views resolve to the same RuntimeG value.
        % G remains the legacy cell schedule; RuntimeG.parameters is ISD.G.
        RuntimeG
        IntBay
        Sampler
        NumSam
        MaxIte
        G                       % optional cell array of prescribed G snapshots
        X
        Y
        FlgCoV
        FlgCvgInt
        CoVTol
        CoVPat
        CoVEvery
        CoVHist
        Verbose
        Estimator               % optional factory: estimator = f(out)
    end
    properties (Access=private)
        StateG = []
    end
    properties (Dependent,Access=private)
        SamplerConfigured
    end

    methods
        function obj = SeMCS(IntBay,Sampler,profile)
            if nargin<1 || nargin>3
                error('SeMIS:Constructor','Provide an IntBay ISD and an optional MCMC sampler.');
            end
            % A bound sampler is already configured. The model-first form
            % retains its existing method-default/profile application.
            samplerFirst=nargin==1 && isobject(IntBay) && isscalar(IntBay) && ...
                isprop(IntBay,'BayStc') && ismethod(IntBay,'SamGen') && ~ismethod(IntBay,'EvlY');
            if samplerFirst
                Sampler=IntBay; IntBay=Sampler.BayStc;
                if isempty(IntBay)
                    error('SeMIS:Constructor','A sampler-first call requires a bound BayStc ISD.');
                end
            elseif nargin<2
                Sampler=[];
            end
            required = {'EvlY','EvlLKF','EvlPDF'};
            for k = 1:numel(required)
                if ~ismethod(IntBay,required{k})
                    error('SeMIS:Interface','IntBay must implement %s.',required{k});
                end
            end
            if ~ismethod(IntBay,'UpdObj') && ~ismethod(IntBay,'UpdObjContext')
                error('SeMIS:Interface','IntBay must implement UpdObj or UpdObjContext.');
            end
            props = {'TF','Ndim','Nfun','G','X','Y','Ncal'};
            for k = 1:numel(props)
                if ~isprop(IntBay,props{k})
                    error('SeMIS:Interface','IntBay must provide %s.',props{k});
                end
            end
            direct=obj.UseDirect(IntBay);
            obj.StateG=SeMISState.Create(IntBay.G,obj.DefaultRun(),Sampler,struct(),IntBay);
            obj.StateG.runtime=struct('iteration',0,'phase','configured','sampler_configured',false);
            space=ISDSpace(IntBay);
            space.Describe(IntBay); % Lift is checked only when MCMC seeds need it.
            if isprop(IntBay,'Nsam') && ~isempty(IntBay.Nsam)
                obj.NumSam = IntBay.Nsam;
            end
            if ismethod(IntBay,'RunDefaults')
                defaults = IntBay.RunDefaults();
                if ~isstruct(defaults) || ~isscalar(defaults)
                    error('SeMIS:RunDefaults','RunDefaults must return a scalar struct.');
                end
                names = fieldnames(defaults);
                for k = 1:numel(names)
                    value=defaults.(names{k});
                    switch names{k}
                        case {'FlgCoV','FlgCvgInt'}
                            valid=islogical(value) && isscalar(value);
                        case 'MaxIte'
                            valid=isnumeric(value) && isreal(value) && isscalar(value) && ...
                                isfinite(value) && value>=1 && value==fix(value);
                        otherwise
                            valid=false;
                    end
                    if ~valid
                        error('SeMIS:RunDefaults','Invalid or unsupported run default: %s.',names{k});
                    end
                    obj.(names{k}) = value;
                end
            end
            if samplerFirst
                obj.SamplerConfigured=true;
            elseif ~direct || ~isempty(Sampler)
                obj.Sampler=obj.ConfigureSampler(IntBay,Sampler);
                obj.SamplerConfigured=true;
            end
            if samplerFirst
                obj.StateG.profile=struct('name','current','evidence', ...
                    struct('provider','bound_sampler','method_class',class(IntBay), ...
                    'sampler_class',class(Sampler)), 'run',struct(),'sampler',struct());
            else
                if nargin<3,profile=[];end
                obj=obj.ConfigureProfile(profile);
            end
        end

        function value=get.RuntimeG(obj), value=obj.StateG; end
        function obj=set.RuntimeG(obj,value)
            value=SeMISState.Normalize(value);
            if isempty(value.model)
                error('SeMIS:RuntimeModel','An engine state requires its full value model.');
            end
            defaults=obj.DefaultRun();names=fieldnames(value.run);
            if ~all(ismember(names,fieldnames(defaults)))
                error('SeMIS:RunState','Unknown runtime run option.');
            end
            for k=1:numel(names),defaults.(names{k})=value.run.(names{k});end
            value.run=defaults;
            obj.ValidateRunOptions(value.run);
            if ~isfield(value.runtime,'sampler_configured'),value.runtime.sampler_configured=false;end
            obj.StateG=value;
        end
        function value=get.IntBay(obj),value=obj.StateG.model;end
        function obj=set.IntBay(obj,value),obj.StateG=SeMISState.WithModel(obj.StateG,value);end
        function value=get.Sampler(obj),value=SeMISState.ApplySampler(obj.StateG);end
        function obj=set.Sampler(obj,value),obj.StateG=SeMISState.CaptureSampler(obj.StateG,value);end
        function value=get.NumSam(obj),value=obj.StateG.run.samples;end
        function obj=set.NumSam(obj,value),obj=obj.SetRunValue('samples',value);end
        function value=get.MaxIte(obj),value=obj.StateG.run.max_iterations;end
        function obj=set.MaxIte(obj,value),obj=obj.SetRunValue('max_iterations',value);end
        function value=get.G(obj),value=obj.StateG.run.schedule;end
        function obj=set.G(obj,value),obj=obj.SetRunValue('schedule',value);end
        function value=get.X(obj),value=obj.StateG.run.x;end
        function obj=set.X(obj,value),obj=obj.SetRunValue('x',value);end
        function value=get.Y(obj),value=obj.StateG.run.y;end
        function obj=set.Y(obj,value),obj=obj.SetRunValue('y',value);end
        function value=get.FlgCoV(obj),value=obj.StateG.run.flag_cov;end
        function obj=set.FlgCoV(obj,value),obj=obj.SetRunValue('flag_cov',value);end
        function value=get.FlgCvgInt(obj),value=obj.StateG.run.flag_isd;end
        function obj=set.FlgCvgInt(obj,value),obj=obj.SetRunValue('flag_isd',value);end
        function value=get.CoVTol(obj),value=obj.StateG.run.cov_tolerance;end
        function obj=set.CoVTol(obj,value),obj=obj.SetRunValue('cov_tolerance',value);end
        function value=get.CoVPat(obj),value=obj.StateG.run.cov_patience;end
        function obj=set.CoVPat(obj,value),obj=obj.SetRunValue('cov_patience',value);end
        function value=get.CoVEvery(obj),value=obj.StateG.run.cov_every;end
        function obj=set.CoVEvery(obj,value),obj=obj.SetRunValue('cov_every',value);end
        function value=get.CoVHist(obj),value=obj.StateG.run.cov_history;end
        function obj=set.CoVHist(obj,value),obj=obj.SetRunValue('cov_history',value);end
        function value=get.Verbose(obj),value=obj.StateG.run.verbose;end
        function obj=set.Verbose(obj,value),obj=obj.SetRunValue('verbose',value);end
        function value=get.Estimator(obj),value=obj.StateG.run.estimator;end
        function obj=set.Estimator(obj,value),obj=obj.SetRunValue('estimator',value);end
        function value=get.SamplerConfigured(obj)
            value=obj.StateG.runtime.sampler_configured;
        end
        function obj=set.SamplerConfigured(obj,value)
            obj.StateG.runtime.sampler_configured=value;
        end

        function obj=ConfigureProfile(obj,name)
            % Profiles belong to the method. Applying one never samples.
            if nargin<2,name=[];end
            if isfield(obj.StateG.runtime,'iteration') && obj.StateG.runtime.iteration>0
                error('SeMIS:ProfileState','Configure a fresh model before starting a different method profile.');
            end
            oldRng=rng;
            try
                spec=SeMISProfile(obj.IntBay,name);
                if ~isequaln(oldRng,rng)
                    error('SeMIS:ProfileRNG','A profile provider must not consume random numbers.');
                end
                candidate=obj;
                if isfield(spec,'model')
                    candidate.IntBay=spec.model;
                    if ismethod(spec.model,'RunDefaults')
                        defaults=spec.model.RunDefaults();
                        if ~isstruct(defaults)||~isscalar(defaults)|| ...
                                ~all(ismember(fieldnames(defaults),{'FlgCoV','FlgCvgInt','MaxIte'}))
                            error('SeMIS:RunDefaults','Invalid profile model run defaults.');
                        end
                        names=fieldnames(defaults);
                        for k=1:numel(names),candidate.(names{k})=defaults.(names{k});end
                    end
                    if ~isempty(candidate.Sampler)
                        candidate.Sampler=candidate.ConfigureSampler(candidate.IntBay,candidate.Sampler);
                    end
                end
                if isfield(spec,'parameters')
                    candidate.StateG=SeMISState.WithParameters(candidate.StateG,spec.parameters);
                end
                if ~all(ismember(fieldnames(spec.run),fieldnames(candidate.DefaultRun())))
                    error('SeMIS:RunState','A profile supplied an unknown run option.');
                end
                candidate.StateG=SeMISState.WithRun(candidate.StateG,spec.run);
                candidate.Sampler=candidate.ApplyProfileSampler(candidate.Sampler,spec.sampler);
                candidate.ValidateRunOptions(candidate.StateG.run);
                if ~isequaln(oldRng,rng)
                    error('SeMIS:ProfileRNG','Applying method and sampler defaults must not consume random numbers.');
                end
                candidate.StateG.profile=struct('name',spec.name,'evidence',spec.evidence, ...
                    'run',spec.run,'sampler',spec.sampler);
                obj=candidate;
            catch cause
                rng(oldRng);rethrow(cause);
            end
        end

        function out = RunIte(obj)
            validateattributes(obj.NumSam,{'numeric'},{'vector','integer','positive','finite'});
            validateattributes(obj.MaxIte,{'numeric'},{'scalar','integer','positive','finite'});
            validateattributes(obj.CoVEvery,{'numeric'},{'scalar','integer','positive','finite'});
            validateattributes(obj.CoVPat,{'numeric'},{'scalar','integer','positive','finite'});
            validateattributes(obj.CoVTol,{'numeric'},{'scalar','positive','finite'});
            intBay = obj.IntBay;
            sampler = obj.Sampler;
            runtimeG=obj.RuntimeG;
            runtimeG.runtime=struct('iteration',0,'phase','initializing', ...
                'sampler_configured',obj.SamplerConfigured,'initial_rng',rng, ...
                'nextPlan',[],'nextSeeds',[],'pending_seed_cost',0);
            samplerConfigured=obj.SamplerConfigured; samplerReset=false;
            maxi = obj.MaxIte;
            vectorBudget = numel(obj.NumSam)>1;
            if vectorBudget
                maxi = min(maxi,numel(obj.NumSam));
                budgets = obj.NumSam(:);
            else
                budgets = repmat(obj.NumSam,maxi,1);
            end
            if ~isempty(obj.G) && ~iscell(obj.G)
                error('SeMIS:Schedule','Prescribed G must be a cell array of ISD snapshots.');
            end
            nPreset = numel(obj.G);
            if nPreset>0, intBay.G = obj.G{1}; end
            if ~obj.UseDirect(intBay)
                if ~samplerConfigured
                    sampler=obj.ConfigureSampler(intBay,sampler); samplerConfigured=true;
                end
                if ismethod(sampler,'ResetState'), sampler=sampler.ResetState; end
                samplerReset=true;
            end
            x = cell(maxi,1); y = cell(maxi,1); g = cell(maxi+1,1);
            xi = cell(maxi,1); yi = cell(maxi,1);
            levels = cell(maxi,1);
            views = cell(maxi,1); descriptors = cell(maxi,1);
            Nsam = zeros(maxi,1); Nsze = zeros(maxi,2); Ncal = zeros(maxi,1);
            MoveShape = zeros(maxi,2); Ngenerated = zeros(maxi,1);
            ChainLengths=cell(maxi,1); SamplerState=cell(maxi,1);
            SeedInfo=cell(maxi,1);
            RuntimeHistory=cell(maxi,1);
            IndependentSamples=false(maxi,1); SamplingMethod=cell(maxi,1);
            Nsze(1,:) = [budgets(1),1];
            g{1} = intBay.G;
            initial = MISInitialState(struct());
            covHist = []; covLevels = [];
            flg_cvg = false;
            reason = 'max_iterations';
            prepareCost = 0;

            for ite = 1:maxi
                runtimeG.runtime.inputPlan=runtimeG.runtime.nextPlan;
                runtimeG.runtime.inputSeeds=runtimeG.runtime.nextSeeds;
                runtimeG.runtime.nextPlan=[];runtimeG.runtime.nextSeeds=[];
                runtimeG.runtime.pending_seed_cost=0;

                % 1. Generate samples in the current ISD's sampling space.
                direct=obj.UseDirect(intBay);
                if direct
                    if ismethod(intBay,'DirectSampleCount')
                        budgets(ite)=intBay.DirectSampleCount(budgets(ite));
                        validateattributes(budgets(ite),{'numeric'},{'scalar','integer','positive','finite'});
                    end
                    sampledG=intBay.G;
                    sampledSize=[intBay.Ndim,intBay.Nfun];
                    [x{ite},y{ite},Ncal(ite),intBay]=intBay.DirectSamples(budgets(ite));
                    obj.CheckDirectSamples(intBay,x{ite},y{ite},sampledG,budgets(ite),sampledSize);
                    if ite==1, initial=obj.DirectInitialState(intBay); end
                    IndependentSamples(ite)=true; SamplingMethod{ite}='direct';
                    Ngenerated(ite)=budgets(ite);
                    SamplerState{ite}=struct('method','direct','independent',true, ...
                        'Nretained',budgets(ite),'Ngenerated',budgets(ite), ...
                        'Ndraws',budgets(ite),'Nproposals',0,'proposalCounts',[], ...
                        'ChainLengths',[]);
                    if ismethod(intBay,'DirectSampleDesign')
                        design=intBay.DirectSampleDesign;
                        if ~isempty(design)
                            design=MISSampleDesign(design,budgets(ite));
                            IndependentSamples(ite)=false;
                            SamplerState{ite}.independent=false;
                            SamplerState{ite}.SampleDesign=design;
                        end
                    end
                elseif ite==1
                    IndependentSamples(ite)=true; SamplingMethod{ite}='initial';
                    sampledG = intBay.G;
                    if ismethod(intBay,'InitialSamples')
                        [x{ite},y{ite},Ncal(ite),intBay,initial] = intBay.InitialSamples(budgets(ite));
                        initial = MISInitialState(struct('initial',initial));
                    else
                        [x{ite},y{ite},Ncal(ite),intBay] = SamGen(intBay,budgets(ite));
                    end
                    if ~isequaln(sampledG,intBay.G)
                        error('SeMIS:InitialSchedule','InitialSamples must not change the requested first ISD.');
                    end
                    Ngenerated(ite) = size(y{ite},2)*size(y{ite},3);
                else
                    SamplingMethod{ite}='mcmc';
                    SeedInfo{ite}=nextSeedInfo;
                    sampler.BayStc = intBay;
                    if ismethod(sampler,'SamGenState')
                        [x{ite},y{ite},Ncal(ite),intBay,sampler,info] = ...
                            sampler.SamGenState(x_sed,y_sed,nextPlan.steps);
                        obj.CheckSamplerInfo(info,y{ite},nextPlan);
                        ChainLengths{ite}=info.ChainLengths;
                        Ngenerated(ite)=info.Ngenerated;
                        SamplerState{ite}=info;
                        if isfield(info,'SeedMode'), SeedInfo{ite}.SeedMode=char(info.SeedMode); end
                    else
                        if ~isscalar(nextPlan.steps)
                            error('SeMIS:RaggedSampler','Unequal chains require a stateful sampler interface.');
                        end
                        [x{ite},y{ite},Ncal(ite),intBay] = ...
                            sampler.SamGen(x_sed,y_sed,nextPlan.steps);
                        Ngenerated(ite)=size(y{ite},2)*size(y{ite},3);
                    end
                    MoveShape(ite,:) = [size(x_sed,2),max(nextPlan.steps)];
                    if strcmp(nextPlan.keep,'last')
                        x{ite} = x{ite}(:,:,end);
                        y{ite} = y{ite}(:,:,end);
                        ChainLengths{ite}=[];
                    end
                end
                if ite==1
                    if initial.estimated && obj.FlgCoV
                        error('SeMIS:InitialNormalizerStopping', ...
                            ['CoV stopping can depend on the estimated initial normalizer. ', ...
                             'Use intrinsic stopping or a fixed level budget to preserve the declared independence.']);
                    end
                    if isfield(initial,'trials') && ~isempty(initial.trials)
                        validateattributes(initial.trials,{'numeric'},{'scalar','integer','positive','finite'});
                        Ngenerated(ite)=initial.trials;
                        if direct
                            SamplerState{ite}.Ngenerated=initial.trials;
                            SamplerState{ite}.Ndraws=initial.trials;
                        end
                    end
                end
                validateattributes(Ncal(ite),{'numeric'},{'scalar','integer','nonnegative','finite'});
                % Seed lifting prepares this sampled level; attribute its
                % model evaluations here, alongside its MCMC evaluations.
                Ncal(ite) = Ncal(ite)+prepareCost;
                prepareCost = 0;
                Nsam(ite) = size(y{ite},2)*size(y{ite},3);
                Nsze(ite,:) = [size(y{ite},2),size(y{ite},3)];
                expectedSamples=budgets(ite);
                if ite>1 && ~direct, expectedSamples=nextPlan.retained; end
                if Nsam(ite)~=expectedSamples
                    error('SeMIS:SampleCount','Sampler returned an unexpected sample count.');
                end
                if size(x{ite},1)~=intBay.Ndim || size(y{ite},1)~=intBay.Nfun || ...
                        ~isequal(size(x{ite},2:3),size(y{ite},2:3))
                    error('SeMIS:SampleLayout','Sample coordinates and raw records have inconsistent dimensions.');
                end
                % 2. Freeze accepted samples, mapped records and runtime state.
                if ismethod(intBay,'AfterSample')
                    sampledG=intBay.G; oldXSize=size(x{ite}); oldYSize=size(y{ite});
                    [x{ite},y{ite},intBay,refreshCost]=intBay.AfterSample(x{ite},y{ite});
                    validateattributes(refreshCost,{'numeric'},{'scalar','integer','nonnegative','finite'});
                    if ~isequaln(sampledG,intBay.G) || ~isequal(oldXSize,size(x{ite})) || ...
                            ~isequal(oldYSize,size(y{ite}))
                        error('SeMIS:RefreshContract','AfterSample must preserve the sampled ISD and record layout.');
                    end
                    Ncal(ite)=Ncal(ite)+refreshCost;
                end
                % Evaluators may hold the last candidate batch; schedules and
                % stopping rules must instead see ALL accepted samples.
                intBay.X = reshape(x{ite},size(x{ite},1),[]);
                intBay.Y = reshape(y{ite}(3:end,:,:),size(y{ite},1)-2,[]);
                space=ISDSpace(intBay);
                projected=~isempty(space.Project);
                if projected
                    sampledG = intBay.G;
                    [xi{ite},yi{ite},intBay,projectCost] = space.ProjectSamples(intBay,x{ite},y{ite});
                    validateattributes(projectCost,{'numeric'},{'scalar','integer','nonnegative','finite'});
                    if ~isequaln(sampledG,intBay.G)
                        error('SeMIS:ProjectionSchedule','Projection must not change the sampled ISD parameters.');
                    end
                    if ~isequal(size(xi{ite},2:3),size(x{ite},2:3)) || ...
                            ~isequal(size(yi{ite},2:3),size(y{ite},2:3)) || size(yi{ite},1)<3
                        error('SeMIS:ProjectionLayout','Projection must preserve the sampled chain layout.');
                    end
                    view = ISDIntegrationView(intBay);
                    if size(xi{ite},1)~=view.Ndim || size(yi{ite},1)~=view.Nfun
                        error('SeMIS:ProjectionLayout','Integration-view dimensions disagree with the projected records.');
                    end
                    Ncal(ite) = Ncal(ite)+projectCost;
                else
                    xi{ite} = x{ite}; yi{ite} = y{ite};
                end
                % Freeze the evaluator AND mapping, including value state
                % outside G, before adaptation can change either one.
                space=ISDSpace(intBay);
                views{ite}=space.IntegrationView(intBay);
                descriptors{ite}=space.Describe(intBay,views{ite});
                if strcmp(space.Kind,'identity') && ~isequaln(x{ite},xi{ite})
                    error('SeMIS:IntegrationContract', ...
                        'An identity adapter must preserve sample coordinates.');
                end
                if ite==1
                    if isfield(initial,'integration') && ...
                            ~isequal(initial.integration,descriptors{ite}.integration)
                        error('SeMIS:IntegrationContract', ...
                            'Initial normalizer metadata uses a different integration space.');
                    end
                    initial.integration=descriptors{ite}.integration;
                    initial.normalizerRelation=descriptors{ite}.normalizerRelation;
                end
                if ite>1 && ~isequal(descriptors{ite}.integration,descriptors{1}.integration)
                    error('SeMIS:IntegrationContract', ...
                        'All levels in one MIS history must use the same integration coordinates and measure.');
                end
                levels{ite}=struct('G',g{ite},'model',intBay,'view',views{ite}, ...
                    'space',space,'descriptor',descriptors{ite}, ...
                    'sampling',struct('x',x{ite},'y',y{ite}), ...
                    'integration',struct('x',xi{ite},'y',yi{ite}), ...
                    'ChainLengths',ChainLengths{ite}, ...
                    'IndependentSamples',IndependentSamples(ite),'Ncal',Ncal(ite), ...
                    'SamplingMethod',SamplingMethod{ite},'SamplerState',SamplerState{ite}, ...
                    'SeedInfo',SeedInfo{ite}, ...
                    'Ngenerated',Ngenerated(ite),'initial',initial);
                RuntimeHistory{ite}=obj.CaptureRuntime(runtimeG,intBay,sampler,ite,'sampled');
                RuntimeHistory{ite}.runtime.requested_budget=budgets(ite);
                RuntimeHistory{ite}.runtime.actual_budget=Nsam(ite);
                RuntimeHistory{ite}.runtime.actual_retained=Nsam(ite);
                RuntimeHistory{ite}.runtime.level_ncal=Ncal(ite);
                levels{ite}.runtimeG=RuntimeHistory{ite};
                dNcal = 0;
                % 3. Form the next ISD without an unused terminal update.
                if ite==maxi
                    % No unused adaptive update after exhausting the budget.
                    g{ite+1} = intBay.G;
                elseif ite+1<=nPreset
                    intBay.G = obj.G{ite+1};
                    g{ite+1} = intBay.G;
                else
                    if ismethod(intBay,'UpdObjContext')
                        [yi{ite},g{ite+1},intBay,dNcal] = intBay.UpdObjContext(levels{ite});
                    else
                        [yi{ite},g{ite+1},intBay,dNcal] = intBay.UpdObj(xi{ite},yi{ite});
                    end
                end
                validateattributes(dNcal,{'numeric'},{'scalar','integer','nonnegative','finite'});
                if ~isequal(size(yi{ite}),size(levels{ite}.integration.y)) || ...
                        ~isequaln(yi{ite}(1:2,:,:),levels{ite}.integration.y(1:2,:,:))
                    error('SeMIS:AdvanceHistory', ...
                        'Adaptation must preserve the sampled density rows and record layout.');
                end
                % Existing adapters may update their accepted record payload
                % while advancing (e.g. aBUS). Identity records stay identical.
                if ~projected, y{ite} = yi{ite}; end
                intBay.G = g{ite+1};
                Ncal(ite) = Ncal(ite)+dNcal;
                levels{ite}.sampling.y=y{ite};
                levels{ite}.integration.y=yi{ite};
                levels{ite}.Ncal=Ncal(ite);
                RuntimeHistory{ite}.runtime.level_ncal=Ncal(ite);
                levels{ite}.runtimeG=RuntimeHistory{ite};
                intBay.X = reshape(x{ite},size(x{ite},1),[]);
                intBay.Y = reshape(y{ite}(3:end,:,:),size(y{ite},1)-2,[]);
                if obj.Verbose, DspRst(ite,Nsam(ite),Nsze(ite,:),Ncal(ite)); end

                runtimeG=obj.CaptureRuntime(runtimeG,intBay,sampler,ite,'advanced');
                context=struct('level',levels{ite},'nextBudget',[], ...
                    'sampler',runtimeG.sampler,'nextG',intBay.G, ...
                    'iteration',ite,'budgetExhausted',ite==maxi,'runtimeG',runtimeG);
                if ite<maxi,context.nextBudget=budgets(ite+1);end
                % 4. Check the method's stopping rule before choosing seeds.
                decision=SeMISMethodPolicy.ShouldStop(intBay,context,'before_seeds');
                if ~isempty(decision) && decision.stop
                    flg_cvg=decision.converged;reason=decision.reason;break
                end

                % Test completion before selecting seeds (a completed
                % schedule must not be reported as an empty-seed failure).
                if isempty(decision) && obj.FlgCvgInt && isprop(intBay,'FlgCvg') && intBay.FlgCvg
                    flg_cvg = true; reason = 'isd_complete'; break
                end
                if obj.FlgCoV && ite>=2 && mod(ite-2,obj.CoVEvery)==0
                    snapshot = struct('x',{x(1:ite)},'y',{y(1:ite)}, ...
                        'g',{g(1:ite+1)},'Nite',ite,'intBay',intBay, ...
                        'Nsam',Nsam(1:ite),'Nsze',Nsze(1:ite,:),'Ncal',Ncal(1:ite));
                    snapshot.integration = struct('x',{xi(1:ite)},'y',{yi(1:ite)}, ...
                        'intBay',views{ite},'views',{views(1:ite)}, ...
                        'descriptors',{descriptors(1:ite)});
                    snapshot.levels=levels(1:ite);
                    snapshot.initial = initial;
                    snapshot.ChainLengths=ChainLengths(1:ite);
                    snapshot.IndependentSamples=IndependentSamples(1:ite);
                    snapshot.SamplingMethod=SamplingMethod(1:ite);
                    snapshot.SamplerState=SamplerState(1:ite);
                    snapshot.SeedInfo=SeedInfo(1:ite);
                    value = obj.EvlCoV(snapshot);
                    covHist(end+1) = value;
                    covLevels(end+1) = ite;
                    if numel(covHist)>=obj.CoVPat+1
                        window = covHist(end-obj.CoVPat:end);
                        if all(isfinite(window)) && all(window>=0) && ...
                                all(diff(window)<=0) && ...
                                abs(window(end)-window(1))/max(abs(window(1)),realmin)<obj.CoVTol
                            flg_cvg = true; reason = 'cov_plateau'; break
                        end
                    end
                end
                if ite==maxi
                    if vectorBudget && maxi==numel(obj.NumSam), reason='sample_budget'; end
                    break
                end
                % Direct samplers draw the next ISD independently; no seeds,
                % chain plan, or conditional lift is needed for that level.
                if obj.UseDirect(intBay)
                    decision=SeMISMethodPolicy.ShouldStop(intBay,context,'after_seeds');
                    if ~isempty(decision) && decision.stop
                        flg_cvg=decision.converged;reason=decision.reason;break
                    end
                    continue
                end
                if ~samplerConfigured
                    sampler=obj.ConfigureSampler(intBay,sampler); samplerConfigured=true;
                end
                if ~samplerReset
                    if ismethod(sampler,'ResetState'), sampler=sampler.ResetState; end
                    samplerReset=true;
                end
                % 5. Select actual parents, then lift them to the next ISD.
                context.sampler=sampler;
                context.sampler.BayStc=intBay;
                [selected,intBay]=SeMISMethodPolicy.SelectSeeds(intBay,context);
                if ~isempty(selected)
                    x_sed=selected.x;y_sed=selected.y;seedCost=selected.cost;
                    nextPlan=selected.plan;nextSeedInfo=selected.seedInfo;
                else
                    plan = [];
                    if ismethod(intBay,'SamplingPlan'), plan = intBay.SamplingPlan(budgets(ite+1)); end
                    if isempty(plan)
                        [x_sed,y_sed,intBay,seedCost,nextSeedInfo] = SedSlt(intBay,xi{ite},yi{ite},budgets(ite+1));
                        Nc = size(x_sed,2);
                        nextPlan = struct('chains',Nc,'steps',1,'keep','all');
                        if Nc>0, nextPlan.steps = budgets(ite+1)/Nc; end
                        nextPlan.retained=budgets(ite+1);
                    else
                        nextPlan = obj.CheckPlan(plan,budgets(ite+1));
                        retainsSeed=isprop(sampler,'SedSav') && sampler.SedSav;
                        if isprop(sampler,'SeedMode'), retainsSeed=retainsSeed || ~strcmp(sampler.SeedMode,'move'); end
                        if strcmp(nextPlan.keep,'last') && retainsSeed
                            error('SeMIS:TransitionCount','Terminal retention requires a sampler that does not count its input seed as a move.');
                        end
                        [x_sed,y_sed,intBay,seedCost,nextSeedInfo] = SedSlt(intBay,xi{ite},yi{ite}, ...
                            budgets(ite+1),nextPlan.chains,nextPlan.selection);
                    end
                end
                prepareCost = seedCost;
                Nc = size(x_sed,2);
                if Nc>0
                    nextSeedInfo.keep=nextPlan.keep;
                    nextSeedInfo.steps=nextPlan.steps;
                    nextSeedInfo.SeedMode='unknown';
                    if isprop(sampler,'SeedMode')
                        nextSeedInfo.SeedMode=char(sampler.SeedMode);
                    elseif isprop(sampler,'SedSav')
                        if sampler.SedSav, nextSeedInfo.SeedMode='retain';
                        else, nextSeedInfo.SeedMode='move'; end
                    end
                end
                runtimeG=obj.CaptureRuntime(runtimeG,intBay,sampler,ite,'seeds_selected');
                runtimeG.runtime.nextPlan=nextPlan;
                runtimeG.runtime.nextSeeds=struct('x',x_sed,'y',y_sed,'info',nextSeedInfo);
                runtimeG.runtime.pending_seed_cost=prepareCost;
                context.runtimeG=runtimeG;context.nextPlan=nextPlan;context.seedInfo=nextSeedInfo;
                context.sampler=runtimeG.sampler;
                decision=SeMISMethodPolicy.ShouldStop(intBay,context,'after_seeds');
                if ~isempty(decision) && decision.stop
                    % A method can stop after lifting seeds. Those actual
                    % evaluations belong to this last retained level.
                    Ncal(ite)=Ncal(ite)+prepareCost;prepareCost=0;
                    levels{ite}.Ncal=Ncal(ite);
                    RuntimeHistory{ite}.runtime.level_ncal=Ncal(ite);
                    levels{ite}.runtimeG=RuntimeHistory{ite};
                    runtimeG.runtime.pending_seed_cost=0;
                    flg_cvg=decision.converged;reason=decision.reason;break
                end
                if Nc==0, reason='empty_support'; break; end
            end
            out = struct('x',{x(1:ite)},'y',{y(1:ite)},'g',{g(1:ite+1)}, ...
                'Nsam',Nsam(1:ite),'Nsze',Nsze(1:ite,:),'Ncal',Ncal(1:ite), ...
                'Nite',ite,'flg_cvg',flg_cvg,'stop_reason',reason, ...
                'intBay',intBay,'CoVHist',covHist,'CoVLevels',covLevels);
            out.initial = initial;
            out.MoveShape = MoveShape(1:ite,:);
            out.Ngenerated = Ngenerated(1:ite);
            out.ChainLengths=ChainLengths(1:ite);
            out.SamplerState=SamplerState(1:ite);
            out.SeedInfo=SeedInfo(1:ite);
            out.IndependentSamples=IndependentSamples(1:ite);
            out.SamplingMethod=SamplingMethod(1:ite);
            out.integration = struct('x',{xi(1:ite)},'y',{yi(1:ite)}, ...
                'intBay',views{ite},'views',{views(1:ite)}, ...
                'descriptors',{descriptors(1:ite)});
            out.levels=levels(1:ite);
            out.runtimeG=RuntimeHistory(1:ite);
            runtimeG=obj.CaptureRuntime(runtimeG,intBay,sampler,ite,'finished');
            runtimeG.runtime.stop_reason=reason;
            runtimeG.runtime.converged=flg_cvg;
            runtimeG.runtime.pending_seed_cost=prepareCost;
            out.finalRuntimeG=runtimeG;
            out.profile=runtimeG.profile;
            out.termination=struct('reason',reason,'converged',flg_cvg, ...
                'budget_exhausted',ismember(reason,{'sample_budget','max_iterations','method_cap'}), ...
                'empty_support',strcmp(reason,'empty_support'));
        end
    end

    methods (Access=private)
        function value=CaptureRuntime(~,base,model,sampler,iteration,phase)
            % Freeze complete value objects and bind the kernel once. Private
            % model/kernel state is retained; density parameters remain opaque.
            value=base;
            value.parameters=model.G; value.model=model; value.sampler=sampler;
            value=SeMISState.Normalize(value);
            value.runtime.iteration=iteration;
            value.runtime.phase=phase;
            value.runtime.rng=rng;
        end
        function obj=SetRunValue(obj,name,value)
            patch=struct();patch.(name)=value;
            obj.StateG=SeMISState.WithRun(obj.StateG,patch);
        end

        function run=DefaultRun(~)
            run=struct('samples',1000,'max_iterations',10000,'schedule',[], ...
                'x',[],'y',[],'flag_cov',true,'flag_isd',false, ...
                'cov_tolerance',.05,'cov_patience',3,'cov_every',1, ...
                'cov_history',[],'verbose',true,'estimator',[]);
        end

        function ValidateRunOptions(~,run)
            try
                validateattributes(run.samples,{'numeric'},{'vector','integer','positive','finite'});
                validateattributes(run.max_iterations,{'numeric'},{'scalar','integer','positive','finite'});
                validateattributes(run.cov_every,{'numeric'},{'scalar','integer','positive','finite'});
                validateattributes(run.cov_patience,{'numeric'},{'scalar','integer','positive','finite'});
                validateattributes(run.cov_tolerance,{'numeric'},{'scalar','positive','finite'});
                for field={'flag_cov','flag_isd','verbose'}
                    value=run.(field{1});
                    if ~((islogical(value)||isnumeric(value)) && isreal(value) && ...
                            isscalar(value) && isfinite(value) && (value==0||value==1))
                        error('SeMIS:RunFlag','%s must be a scalar logical or 0/1.',field{1});
                    end
                end
                if ~isempty(run.schedule) && ~iscell(run.schedule)
                    error('SeMIS:RunSchedule','A prescribed density schedule must be a cell array.');
                end
                if ~isempty(run.estimator) && ~isa(run.estimator,'function_handle')
                    error('SeMIS:RunEstimator','Estimator must be a function handle or empty.');
                end
            catch cause
                error('SeMIS:RunState','Invalid runtime run options: %s',cause.message);
            end
        end

        function sampler=ApplyProfileSampler(obj,sampler,patch)
            if isempty(fieldnames(patch)),return;end
            if isempty(sampler)
                error('SeMIS:ProfileSampler','A sampler profile needs an explicitly supplied kernel.');
            end
            sampler=obj.ConfigureSampler(obj.IntBay,sampler,patch);
        end

        function direct = UseDirect(~,intBay)
            direct=ismethod(intBay,'DirectSamples');
            if ismethod(intBay,'DirectSampling')
                direct=intBay.DirectSampling();
                if ~islogical(direct) || ~isscalar(direct)
                    error('SeMIS:DirectSampling','DirectSampling must return a logical scalar.');
                end
                if direct && ~ismethod(intBay,'DirectSamples')
                    error('SeMIS:DirectSampling','An enabled direct ISD must implement DirectSamples.');
                end
            end
        end

        function initial = DirectInitialState(~,intBay)
            if ismethod(intBay,'DirectInitialState')
                initial=MISInitialState(struct('initial',intBay.DirectInitialState()));
                return
            end
            value=[];
            if ismethod(intBay,'KnownLogNormalizer'), value=intBay.KnownLogNormalizer(); end
            if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
                error('SeMIS:DirectInitialNormalizer', ...
                    ['A direct first ISD needs a finite KnownLogNormalizer or ', ...
                     'explicit DirectInitialState normalizer metadata.']);
            end
            initial=MISInitialState(struct()); initial.logC=value;
            initial.estimator='known_direct_isd';
        end

        function CheckDirectSamples(~,intBay,u,y,sampledG,N,sampledSize)
            if ~isequaln(sampledG,intBay.G)
                error('SeMIS:DirectSchedule','DirectSamples must preserve the sampled G.');
            end
            if ~isequal([intBay.Ndim,intBay.Nfun],sampledSize)
                error('SeMIS:DirectLayout','DirectSamples must preserve sampling dimensions.');
            end
            if ~isnumeric(u) || ~isreal(u) || ~isequal(size(u),[intBay.Ndim,N]) || ...
                    any(~isfinite(u),'all') || ~isnumeric(y) || ~isreal(y) || ...
                    ~isequal(size(y),[intBay.Nfun,N])
                error('SeMIS:DirectLayout','DirectSamples must return N independent columns of coordinates and complete raw records.');
            end
            probe=intBay; probe.X=u; probe.Y=y(3:end,:);
            fresh=[probe.EvlLKF;probe.EvlPDF];
            if ~isnumeric(fresh) || ~isreal(fresh) || ~isequal(size(fresh),[2,N]) || ...
                    any(~isfinite(fresh),'all') || ...
                    any(~isfinite(y(1:2,:)),'all') || ...
                    any(abs(fresh-y(1:2,:))>1e-10*(1+abs(fresh)),'all')
                error('SeMIS:DirectDensity','DirectSamples raw density must match its current ISD.');
            end
        end

        function Sampler = ConfigureSampler(~,IntBay,Sampler,defaults)
            if ~isprop(Sampler,'BayStc') || ~ismethod(Sampler,'SamGen')
                error('SeMIS:Interface','A non-direct ISD requires a sampler with BayStc and SamGen.');
            end
            if nargin>=4 || ismethod(IntBay,'SamplerDefaults')
                if nargin<4,defaults = IntBay.SamplerDefaults();end
                if ~isstruct(defaults) || ~isscalar(defaults)
                    error('SeMIS:SamplerDefaults','SamplerDefaults must return a scalar struct.');
                end
                names = fieldnames(defaults);
                for k = 1:numel(names)
                    name = names{k};
                    if ~ismember(name,{'InitialCorrelation','InitialScale','AdaptChains','AdaptProposals','SeedMode', ...
                            'ScaleMode','CarryScale','AdaptStatistic','AdaptFraction','AdaptPartial', ...
                            'AdaptParameter','AdaptExponent','MaxCorrelation', ...
                            'MixingMode','MixingProbeFactory','MixingGain','MixingMinProposals', ...
                            'MixingMinJumps','MixingMaxRatio'})
                        error('SeMIS:SamplerDefaults','Unsupported sampler policy: %s.',name);
                    end
                    value = defaults.(name);
                    switch name
                        case {'SeedMode','ScaleMode','AdaptStatistic','MixingMode','AdaptParameter'}
                            choices=struct('SeedMode',{{'move','retain','reevaluate'}}, ...
                                'ScaleMode',{{'unit','seeds'}},'AdaptStatistic',{{'realized','probability'}}, ...
                                'AdaptParameter',{{'lambda','rho'}}, ...
                                'MixingMode',{{'off','diagnostic','score'}});
                            valid=ischar(value) && isrow(value) && ismember(value,choices.(name));
                        case 'MixingProbeFactory'
                            valid=isempty(value) || isa(value,'function_handle');
                        case {'CarryScale','AdaptPartial'}
                            valid=islogical(value) && isscalar(value);
                        otherwise
                            valid=isnumeric(value) && isreal(value) && ...
                                (isempty(value) || (isscalar(value) && isfinite(value)));
                            if valid && ~isempty(value)
                                if strcmp(name,'InitialCorrelation')
                                    valid=value>=0 && value<=1;
                                elseif ismember(name,{'AdaptFraction','InitialScale'})
                                    valid=value>0;
                                elseif strcmp(name,'AdaptExponent')
                                    valid=value>0 && value<=1;
                                elseif strcmp(name,'MaxCorrelation')
                                    valid=value>0 && value<1;
                                elseif strcmp(name,'MixingGain')
                                    valid=value>=0;
                                elseif strcmp(name,'MixingMaxRatio')
                                    valid=value>=1;
                                else
                                    valid=value>=1 && value==fix(value);
                                end
                            end
                            if ismember(name,{'MixingGain','MixingMinProposals','MixingMinJumps','MixingMaxRatio', ...
                                    'AdaptExponent','MaxCorrelation'})
                                valid=valid && ~isempty(value);
                            end
                    end
                    if ~valid
                        error('SeMIS:SamplerDefaults','Invalid value for sampler policy %s.',name);
                    end
                end
                % Kernels may interpret valid ISD policies through their
                % own capability hook. For example, an elliptical slice
                % kernel has no aCS proposal scale or acceptance target.
                if ismethod(Sampler,'ConfigureISDPolicy')
                    configured = Sampler.ConfigureISDPolicy(defaults);
                    if ~isscalar(configured) || ~strcmp(class(configured),class(Sampler)) || ...
                            ~isprop(configured,'BayStc') || ~ismethod(configured,'SamGen')
                        error('SeMIS:SamplerDefaults', ...
                            'ConfigureISDPolicy must return a scalar sampler of the same class.');
                    end
                    Sampler = configured;
                    return
                end
                samplerMeta = metaclass(Sampler);
                propertiesMeta = samplerMeta.PropertyList;
                for k = 1:numel(names)
                    name = names{k};
                    property = propertiesMeta(strcmp({propertiesMeta.Name},name));
                    if isempty(property) || property.Constant || ...
                            ~ischar(property.SetAccess) || ~strcmp(property.SetAccess,'public')
                        error('SeMIS:SamplerDefaults', ...
                            'Sampler must expose a publicly writable %s property.',name);
                    end
                end
                % Validate the whole policy before changing the sampler.
                for k = 1:numel(names)
                    name = names{k};
                    try
                        Sampler.(name) = defaults.(name);
                    catch cause
                        error('SeMIS:SamplerDefaults','Cannot set sampler policy %s: %s',name,cause.message);
                    end
                end
            end
        end

        function CheckSamplerInfo(~,info,y,plan)
            fields={'ChainLengths','Nretained','Ngenerated','Nproposals','proposalCounts'};
            if ~isstruct(info) || ~isscalar(info) || ~all(isfield(info,fields))
                error('SeMIS:SamplerMetadata','SamGenState must report complete sample and chain metadata.');
            end
            for name={'Nretained','Ngenerated','Nproposals'}
                validateattributes(info.(name{1}),{'numeric'},{'scalar','integer','nonnegative','finite'});
            end
            validateattributes(info.proposalCounts,{'numeric'},{'vector','integer','nonnegative','finite'});
            generatedPositions=plan.chains*plan.steps;
            if ~isscalar(plan.steps), generatedPositions=sum(plan.steps); end
            if info.Nretained~=size(y,2)*size(y,3) || info.Nretained~=generatedPositions || ...
                    numel(info.proposalCounts)~=plan.chains || sum(info.proposalCounts)~=info.Nproposals || ...
                    info.Ngenerated~=info.Nproposals
                error('SeMIS:SamplerMetadata','Sample counts and sampler metadata disagree.');
            end
            if isscalar(plan.steps)
                if ~isempty(info.ChainLengths) || ~isequal(size(y,2:3),[plan.chains,plan.steps])
                    error('SeMIS:SamplerMetadata','A rectangular plan needs a rectangular output layout.');
                end
            elseif ~isequal(info.ChainLengths(:),plan.steps(:)) || size(y,3)~=1
                error('SeMIS:SamplerMetadata','Packed output must report the actual ordered chain lengths.');
            end
        end

        function plan = CheckPlan(~,plan,budget)
            if ~isstruct(plan) || ~isscalar(plan) || ~all(isfield(plan,{'chains','steps','keep'}))
                error('SeMIS:SamplingPlan','SamplingPlan must return chains, steps and keep.');
            end
            validateattributes(plan.chains,{'numeric'},{'scalar','integer','positive','finite'});
            validateattributes(plan.steps,{'numeric'},{'vector','integer','positive','finite'});
            if ~isscalar(plan.steps) && numel(plan.steps)~=plan.chains
                error('SeMIS:SamplingPlan','A step vector needs one length per chain.');
            end
            plan.keep = validatestring(plan.keep,{'all','last'});
            if ~isfield(plan,'budgetMode'), plan.budgetMode='exact'; end
            plan.budgetMode=validatestring(plan.budgetMode,{'exact','nominal'});
            if ~isfield(plan,'selection'), plan.selection='weighted'; end
            plan.selection=validatestring(plan.selection,{'weighted','survivors'});
            retained = plan.chains;
            if strcmp(plan.keep,'all')
                if isscalar(plan.steps), retained=retained*plan.steps; else, retained=sum(plan.steps); end
            elseif ~isscalar(plan.steps)
                error('SeMIS:SamplingPlan','Terminal retention currently requires equal chain lengths.');
            end
            if retained~=budget && strcmp(plan.budgetMode,'exact')
                error('SeMIS:SamplingPlan','The plan must retain exactly the requested sample budget.');
            end
            plan.retained=retained;
        end

        function value = EvlCoV(obj,out)
            if ~isempty(obj.Estimator)
                mis = obj.Estimator(out);
            elseif ismethod(out.intBay.TF,'EvlLSF') || ...
                    (isstruct(out.intBay.TF) && isfield(out.intBay.TF,'EvlLSF') && ...
                     isa(out.intBay.TF.EvlLSF,'function_handle'))
                mis = MIS_Rel(out);
            else
                mis = MIS_Bay(out);
            end
            [~,value] = mis.estimateCov;
            value = max(value(:));
            if ~isreal(value) || isempty(value) || any(~isfinite(value))
                value = inf;
            end
        end
    end
end
