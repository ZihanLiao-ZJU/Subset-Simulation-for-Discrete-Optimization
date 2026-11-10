function initial = MISInitialState(in)
% Validate the normalizer shared by every proposal in a MIS history.
% A missing metadata record preserves the legacy normalized-reference API.
% A stochastic initial normalizer is supported only when independent of the
% retained conditional samples and subsequent adaptation/stopping decisions.
if ~isfield(in,'initial')
    initial = struct('logC',0,'cov',0,'independent',true,'estimated',false, ...
        'estimator','known_normalized_reference','trials',[],'successes',[]);
    return
end
initial = in.initial;
assert(isstruct(initial) && isscalar(initial) && ...
    isfield(initial,'logC') && isfield(initial,'cov'), ...
    'SeMIS:InitialNormalizer','Initial metadata requires scalar logC and cov fields.');
assert(isnumeric(initial.logC) && isscalar(initial.logC) && ...
    isreal(initial.logC) && isfinite(initial.logC), ...
    'SeMIS:InitialNormalizer','The initial log normalizer must be a finite real scalar.');
assert(isnumeric(initial.cov) && isscalar(initial.cov) && ...
    isreal(initial.cov) && isfinite(initial.cov) && initial.cov>=0, ...
    'SeMIS:InitialNormalizer','The initial normalizer CoV must be finite and nonnegative.');
if isfield(initial,'estimated')
    initial.estimated = logicalScalar(initial.estimated,'estimated');
else
    initial.estimated = initial.cov>0;
end
assert(initial.estimated || initial.cov==0, ...
    'SeMIS:InitialNormalizer','A positive initial normalizer CoV requires estimated=true.');
if isfield(initial,'independent')
    initial.independent = logicalScalar(initial.independent,'independent');
elseif initial.estimated
    error('SeMIS:InitialDependence', ...
        'An estimated initial normalizer requires explicit independent=true metadata, even when its plug-in CoV is zero.');
else
    initial.independent = true;
end
if initial.estimated && ~initial.independent
    error('SeMIS:InitialDependence', ...
        'Dependence between the initial normalizer and the MIS history is unsupported.');
end
% Retain provenance and any additional caller fields without reinterpretation.
if ~isfield(initial,'estimator'), initial.estimator = 'unspecified'; end
if ~isfield(initial,'trials'), initial.trials = []; end
if ~isfield(initial,'successes'), initial.successes = []; end
end

function value = logicalScalar(value,name)
assert((islogical(value) || isnumeric(value)) && isscalar(value) && ...
    isreal(value) && isfinite(value) && (value==0 || value==1), ...
    'SeMIS:InitialNormalizer','Initial %s must be a logical scalar.',name);
value = logical(value);
end
