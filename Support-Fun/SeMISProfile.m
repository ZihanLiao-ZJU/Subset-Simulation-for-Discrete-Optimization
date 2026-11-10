function spec = SeMISProfile(model,name)
% Read and validate one method-owned profile; do not apply its option patches.
% No method-family registry, sampler construction, target evaluation or RNG.
% When both model and parameters are provided, the caller applies the model
% first, then the explicit parameters patch through SeMISState.WithParameters.
if nargin<2,name=[];end
emptyName=isempty(name) || (isstring(name) && isscalar(name) && strlength(name)==0);
if ~emptyName
    if isstring(name) && isscalar(name),name=char(name);end
    if ~ischar(name) || ~isrow(name)
        error('SeMIS:MethodProfile','Profile name must be a nonempty character row or scalar string.');
    end
end
if ~isobject(model) || ~isscalar(model) || isa(model,'handle') || ~isprop(model,'G')
    error('SeMIS:MethodProfile','Profiles require a scalar value model with a legacy G property.');
end
if ismethod(model,'MethodProfileSpec')
    if emptyName,spec=model.MethodProfileSpec();else,spec=model.MethodProfileSpec(name);end
else
    if ~emptyName
        error('SeMIS:ProfileProvider','An explicit profile requires model.MethodProfileSpec(name).');
    end
    protocol=[];
    if isprop(model,'Protocol'),protocol=model.Protocol;end
    profileName='current';
    if ischar(protocol) && isrow(protocol) && ~isempty(protocol),profileName=protocol;
    elseif isstring(protocol) && isscalar(protocol) && strlength(protocol)>0,profileName=char(protocol);end
    spec=struct('name',profileName,'run',struct(),'sampler',struct(), ...
        'evidence',struct('provider','legacy_current_settings','method_class',class(model), ...
        'protocol',protocol,'metadata_only',true));
end
if ~isstruct(spec) || ~isscalar(spec) || ~isfield(spec,'name')
    error('SeMIS:MethodProfile','MethodProfileSpec must return a scalar struct with a canonical name.');
end
allowed={'name','run','sampler','evidence','parameters','model'};
if ~all(ismember(fieldnames(spec),allowed))
    error('SeMIS:MethodProfile','Unknown profile-spec field; put provenance metadata in evidence.');
end
if isstring(spec.name) && isscalar(spec.name),spec.name=char(spec.name);end
if ~ischar(spec.name) || ~isrow(spec.name) || isempty(spec.name)
    error('SeMIS:MethodProfile','Profile spec must supply a nonempty canonical name.');
end
for field={'run','sampler','evidence'}
    key=field{1};
    if ~isfield(spec,key),spec.(key)=struct();end
    if ~isstruct(spec.(key)) || ~isscalar(spec.(key))
        error('SeMIS:MethodProfile','Profile %s patch must be a scalar struct.',key);
    end
end
if isfield(spec,'model')
    candidate=spec.model;
    if ~isobject(candidate) || ~isscalar(candidate) || isa(candidate,'handle') || ...
            ~strcmp(class(candidate),class(model)) || ~isprop(candidate,'G')
        error('SeMIS:ProfileModel','A profile model must be a scalar value object of the same class.');
    end
    if isprop(model,'TF') && (~isprop(candidate,'TF') || ~isequaln(candidate.TF,model.TF))
        error('SeMIS:ProfileTarget','A profile must preserve the target TF; it may not substitute a new problem.');
    end
end
end
