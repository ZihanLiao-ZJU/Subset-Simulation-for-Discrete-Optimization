function x = U2X(intBay,u)
% Dist is authoritative when nonempty; Pdis is the legacy model interface.
% Never hide a malformed Dist by silently switching to a different prior.
target = intBay.TF;
if hasMember(target,'Dist') && ~isempty(target.Dist)
    transform = target.Dist;
elseif hasMember(target,'Pdis') && ~isempty(target.Pdis)
    transform = target.Pdis;
else
    error('SeMIS:PriorTransform','The model needs a nonempty Dist or Pdis transform.');
end
callable = isscalar(transform) && ...
    ((isobject(transform) && ismethod(transform,'U2X')) || ...
    (hasMember(transform,'U2X') && isa(transform.U2X,'function_handle')));
if ~callable
    error('SeMIS:PriorTransform','The selected prior transform must expose callable U2X.');
end
x = transform.U2X(u);
end

function found = hasMember(value,name)
found = (isstruct(value) && isscalar(value) && isfield(value,name)) || ...
    (isobject(value) && isscalar(value) && isprop(value,name));
end
