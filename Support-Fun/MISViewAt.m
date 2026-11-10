function view = MISViewAt(evaluators,g,k)
% A modern history freezes the complete value evaluator at every level.
% A single legacy evaluator still takes its historical parameters from G.
if iscell(evaluators)
    assert(k<=numel(evaluators) && ~isempty(evaluators{k}), ...
        'SeMIS:IntegrationHistory','A frozen integration view is required for every sampled level.');
    view=evaluators{k};
    assert(isequaln(view.G,g{k}),'SeMIS:IntegrationHistory', ...
        'The frozen integration view and recorded proposal parameters disagree.');
else
    view=evaluators;
    view.G=g{k};
end
assert(~isa(view,'handle'),'SeMIS:MutableHistory', ...
    'Historical integration evaluators must be value objects.');
end
