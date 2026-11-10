function space = ISDSpace(model)
% Resolve an explicit value adapter, legacy projection hooks, or identity.
assert(~isa(model,'handle'),'SeMIS:MutableHistory', ...
    'ISDs must be value objects; mutable handle state cannot freeze MIS history.');
if isprop(model,'SpaceAdapter') && ~isempty(model.SpaceAdapter)
    space=model.SpaceAdapter;
elseif ismethod(model,'GetSpaceAdapter')
    space=model.GetSpaceAdapter();
else
    space=ISDSpaceAdapter;
    hooks=[ismethod(model,'IntegrationView'),ismethod(model,'ProjectSamples')];
    if any(hooks)
        assert(all(hooks),'SeMIS:ProjectionInterface', ...
            'Legacy projection requires IntegrationView and ProjectSamples together.');
        space.Kind='legacy';
        space.IntegrationCoordinates='legacy-integration';
        space.View=@(m)m.IntegrationView();
        space.Project=@legacyProject;
        if ismethod(model,'LiftSamples'), space.Lift=@legacyLift; end
    elseif ismethod(model,'LiftSamples')
        error('SeMIS:ProjectionInterface','LiftSamples requires an integration projection.');
    end
end
assert(isa(space,'ISDSpaceAdapter') && isscalar(space), ...
    'SeMIS:SpaceAdapter','SpaceAdapter/GetSpaceAdapter must supply an ISDSpaceAdapter value.');
% An explicit adapter used only for output must not silently re-enable
% row-based diagnostics for an unrelated custom target or payload.
customPayload=(isprop(model,'LogTarget') && ~isempty(model.LogTarget)) || ...
    (isprop(model,'EvaluatePayload') && ~isempty(model.EvaluatePayload));
if customPayload && isempty(space.LikelihoodDiagnostics)
    space.LikelihoodDiagnostics=false;
end
end

function [x,y,m,cost]=legacyProject(m,x,y)
[x,y,m,cost]=m.ProjectSamples(x,y);
end
function [x,y,m,cost]=legacyLift(m,x,y)
[x,y,m,cost]=m.LiftSamples(x,y);
end
