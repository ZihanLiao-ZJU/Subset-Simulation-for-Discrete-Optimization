classdef SeMISState
    % Value-state envelope around an unchanged legacy density payload.
    % parameters is exactly model.G: numeric/struct/cell payloads are opaque.
    % Full value objects retain private/protected model and kernel evolution.
    % This is not a resume API. RNG, next seeds/plan and history position must
    % be recorded explicitly by the caller; no method here reads/restores RNG.
    methods (Static)
        function G = Create(parameters,run,sampler,profile,model)
            if nargin<2,run=struct();end
            if nargin<3,sampler=[];end
            if nargin<4,profile=struct();end
            if nargin<5,model=[];end
            G=struct('schema','semis.runtime.v1','parameters',[], ...
                'run',run,'sampler',[],'model',[],'profile',profile,'runtime',struct());
            % Assignment after construction avoids expanding a cell payload
            % into a struct array through MATLAB's struct constructor.
            G.parameters=parameters;G.sampler=sampler;G.model=model;
            G=SeMISState.Normalize(G);
        end

        function G = Normalize(G)
            SeMISState.CheckEnvelope(G);
            if ~isfield(G,'model'),G.model=[];end
            if ~isfield(G,'runtime'),G.runtime=struct();end
            SeMISState.CheckParts(G);
            if ~isempty(G.model)
                model=G.model;
                model.G=G.parameters;
                G.model=model;
            end
            if ~isempty(G.sampler) && ~isempty(G.model) && isprop(G.sampler,'BayStc')
                sampler=G.sampler;
                sampler.BayStc=G.model;
                G.sampler=sampler;
            end
            G=SeMISState.Validate(G);
        end

        function G = Validate(G)
            SeMISState.CheckEnvelope(G);
            if ~isfield(G,'model') || ~isfield(G,'runtime')
                error('SeMIS:RuntimeState','Runtime state requires model and runtime fields; use Normalize for defaults.');
            end
            SeMISState.CheckParts(G);
            if ~isempty(G.model) && ~isequaln(G.model.G,G.parameters)
                error('SeMIS:RuntimeStateSync','Runtime parameters and the legacy model.G view disagree.');
            end
        end

        function G = WithParameters(G,parameters)
            G=SeMISState.Validate(G);
            G.parameters=parameters;
            G=SeMISState.Normalize(G);
        end

        function G = WithModel(G,model)
            G=SeMISState.Validate(G);
            SeMISState.CheckModel(model);
            if isempty(model)
                error('SeMIS:RuntimeModel','WithModel requires a scalar value model.');
            end
            G.model=model;G.parameters=model.G;
            G=SeMISState.Normalize(G);
        end

        function G = WithRun(G,patch)
            G=SeMISState.Validate(G);
            SeMISState.CheckStruct(patch,'run patch');
            names=fieldnames(patch);
            for k=1:numel(names),G.run.(names{k})=patch.(names{k});end
            G=SeMISState.Validate(G);
        end

        function G = CaptureSampler(G,sampler)
            % Call WithModel first when a sampling call returned a new model.
            % The kernel's BayStc is a bound view, not another authority.
            G=SeMISState.Validate(G);
            SeMISState.CheckSampler(sampler);
            G.sampler=sampler;
            G=SeMISState.Normalize(G);
        end

        function sampler = ApplySampler(G)
            G=SeMISState.Validate(G);
            sampler=G.sampler;
            if ~isempty(sampler) && ~isempty(G.model) && isprop(sampler,'BayStc')
                sampler.BayStc=G.model;
            end
        end
    end

    methods (Static, Access=private)
        function CheckEnvelope(G)
            if ~isstruct(G) || ~isscalar(G) || ~isfield(G,'schema') || ...
                    ~(ischar(G.schema) && isrow(G.schema) && strcmp(G.schema,'semis.runtime.v1'))
                error('SeMIS:RuntimeState','Expected a scalar semis.runtime.v1 state envelope.');
            end
            required={'parameters','run','sampler','profile'};
            if ~all(isfield(G,required))
                error('SeMIS:RuntimeState','Runtime state is missing a required envelope field.');
            end
        end

        function CheckParts(G)
            SeMISState.CheckStruct(G.run,'run');
            SeMISState.CheckStruct(G.profile,'profile');
            SeMISState.CheckStruct(G.runtime,'runtime');
            SeMISState.CheckModel(G.model);
            SeMISState.CheckSampler(G.sampler);
        end

        function CheckStruct(value,name)
            if ~isstruct(value) || ~isscalar(value)
                error('SeMIS:RuntimeState','%s must be a scalar struct.',name);
            end
        end

        function CheckModel(model)
            if isnumeric(model) && isempty(model),return;end
            if ~isobject(model) || ~isscalar(model) || isa(model,'handle') || ~isprop(model,'G')
                error('SeMIS:RuntimeModel','Model must be [] or a scalar value object with a legacy G property.');
            end
        end

        function CheckSampler(sampler)
            if isnumeric(sampler) && isempty(sampler),return;end
            if ~isobject(sampler) || ~isscalar(sampler) || isa(sampler,'handle')
                error('SeMIS:RuntimeSampler','Sampler must be [] or a scalar value kernel object.');
            end
        end
    end
end
