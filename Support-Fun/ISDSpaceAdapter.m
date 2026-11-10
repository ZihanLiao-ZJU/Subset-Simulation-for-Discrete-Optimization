classdef ISDSpaceAdapter
    % Value adapter separating sampling, integration and posterior output.
    % Callbacks must use value state, not mutable external/handle captures.
    % A projection preserves total unnormalized mass, including Jacobians.
    properties
        Kind = 'identity' % identity | bijective | marginal | legacy
        SamplingCoordinates = 'u'
        IntegrationCoordinates = 'u'
        SamplingMeasure = 'lebesgue'
        IntegrationMeasure = 'lebesgue'
        OutputCoordinates = 'physical'
        OutputMeasure = 'pushforward'
        OutputDimension = []
        NormalizerRelation = 'preserved'
        LikelihoodDiagnostics = [] % [] enables legacy diagnostics only for identity
        View = []       % @(model) integration density evaluator (value object)
        Project = []    % @(model,x,y) -> [xi,yi,model,cost]
        Lift = []       % @(model,xi,yi) -> [x,y,model,cost], NEXT proposal conditional
        Posterior = []  % @(sampledModel,xi,yi) -> [x,y], TARGET conditional
        Output = []     % @(sampledModel,xi,yi) -> [x,y], user output mapping
    end
    methods
        function view = IntegrationView(obj,model)
            if isempty(obj.View)
                assert(strcmp(obj.Kind,'identity'),'SeMIS:ProjectionInterface', ...
                    'A nonidentity space requires a View callback.');
                view=model;
            else
                view=obj.View(model);
            end
            assert(~isa(view,'handle'),'SeMIS:MutableHistory', ...
                'Integration evaluators must be value objects so historical densities are frozen.');
        end

        function [x,y,model,cost] = ProjectSamples(obj,model,x,y)
            cost=0;
            if ~isempty(obj.Project)
                [x,y,model,cost]=obj.Project(model,x,y);
            elseif ~strcmp(obj.Kind,'identity')
                error('SeMIS:ProjectionInterface','A nonidentity space requires a Project callback.');
            end
        end

        function [x,y,model,cost] = LiftSamples(obj,model,x,y)
            cost=0;
            if ~isempty(obj.Lift)
                [x,y,model,cost]=obj.Lift(model,x,y);
            elseif strcmp(obj.Kind,'identity') && isempty(obj.Project)
                y=UpdY(model,x,y);
            else
                error('SeMIS:MissingProposalLift', ...
                    'MCMC seeds require a proposal-conditional Lift when coordinates or projected payload may change; direct sampling does not.');
            end
        end

        function [x,y] = PosteriorSamples(obj,model,x,y)
            if ~isempty(obj.Posterior)
                [x,y]=obj.Posterior(model,x,y);
            elseif ~strcmp(obj.Kind,'identity')
                error('SeMIS:MissingPosteriorMap', ...
                    'Full posterior recovery requires a target-conditional Posterior callback, not proposal Lift.');
            end
        end

        function [x,y] = OutputSamples(obj,model,x,y)
            if ~isempty(obj.Output)
                [x,y]=obj.Output(model,x,y);
            elseif strcmp(obj.Kind,'identity') && size(x,1)==model.TF.Ndim && ...
                    strcmp(obj.IntegrationCoordinates,'u') && ...
                    strcmp(obj.IntegrationMeasure,'lebesgue') && ...
                    strcmp(obj.OutputCoordinates,'physical') && ...
                    strcmp(obj.OutputMeasure,'pushforward')
                x=U2X(model,x);
            else
                error('SeMIS:MissingOutputMap', ...
                    'This space needs an explicit Output callback to return physical parameters.');
            end
        end

        function descriptor = Describe(obj,model,view)
            if nargin<3, view=obj.IntegrationView(model); end
            assert(any(strcmp(obj.Kind,{'identity','bijective','marginal','legacy'})), ...
                'SeMIS:SpaceDescriptor','Kind must be identity, bijective, marginal or legacy.');
            callbacks={'View','Project','Lift','Posterior','Output'};
            for k=1:numel(callbacks)
                value=obj.(callbacks{k});
                assert(isempty(value) || isa(value,'function_handle'), ...
                    'SeMIS:SpaceCallback','Space callbacks must be function handles or empty.');
            end
            if ~strcmp(obj.Kind,'identity')
                assert(~isempty(obj.View) && ~isempty(obj.Project), ...
                    'SeMIS:ProjectionInterface','A nonidentity space requires View and Project callbacks.');
            end
            assert(strcmp(obj.NormalizerRelation,'preserved'), ...
                'SeMIS:SpaceNormalizer', ...
                'Space mappings must preserve the unnormalized integral; include Jacobians or marginalization factors.');
            fields={'SamplingCoordinates','IntegrationCoordinates','SamplingMeasure', ...
                'IntegrationMeasure','OutputCoordinates','OutputMeasure'};
            for k=1:numel(fields)
                value=obj.(fields{k});
                assert((ischar(value) && isrow(value) && ~isempty(value)) || ...
                    (isstring(value) && isscalar(value) && strlength(value)>0), ...
                    'SeMIS:SpaceDescriptor','Space names and measures must be nonempty scalar text.');
            end
            outputDim=obj.OutputDimension;
            if isempty(outputDim), outputDim=model.TF.Ndim; end
            dims=[model.Ndim,view.Ndim,outputDim];
            validateattributes(dims,{'numeric'},{'vector','integer','positive','finite'});
            diagnostics=obj.LikelihoodDiagnostics;
            if isempty(diagnostics)
                % A new evaluator/projector may reinterpret cached payload
                % even when the coordinates themselves are unchanged.
                diagnostics=strcmp(obj.Kind,'identity') && isempty(obj.View) && isempty(obj.Project);
            end
            assert(islogical(diagnostics) && isscalar(diagnostics), ...
                'SeMIS:SpaceDescriptor','LikelihoodDiagnostics must be a logical scalar.');
            descriptor=struct('kind',char(obj.Kind), ...
                'sampling',struct('coordinates',char(obj.SamplingCoordinates), ...
                    'dimension',model.Ndim,'measure',char(obj.SamplingMeasure)), ...
                'integration',struct('coordinates',char(obj.IntegrationCoordinates), ...
                    'dimension',view.Ndim,'measure',char(obj.IntegrationMeasure)), ...
                'output',struct('coordinates',char(obj.OutputCoordinates), ...
                    'dimension',outputDim,'measure',char(obj.OutputMeasure)), ...
                'normalizerRelation',char(obj.NormalizerRelation), ...
                'likelihoodDiagnostics',diagnostics);
            if strcmp(obj.Kind,'identity')
                assert(isequal(descriptor.sampling,descriptor.integration), ...
                    'SeMIS:SpaceDescriptor','Identity sampling and integration spaces must agree.');
            end
        end
    end
end
