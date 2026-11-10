function logMass = MisParticleWeights(w_mis,w_target)
% A level contributes mean(W_i), so each pooled particle has mass W_ik/N_i.
% Omitting /N_i overrepresents levels with more samples in posterior draws.
assert(numel(w_mis)==numel(w_target),'SeMIS:MISSize','Weight histories differ.');
parts = cell(numel(w_mis),1);
for k = 1:numel(parts)
    assert(isequal(size(w_mis{k}),size(w_target{k})), ...
        'SeMIS:MISSize','Weights must preserve [number of chains, chain length].');
    parts{k} = reshape(w_mis{k}+w_target{k},1,[]) - log(numel(w_mis{k}));
end
logMass = [parts{:}];
normalizer = logsum(logMass,'all');
assert(isreal(normalizer) && isfinite(normalizer), ...
    'SeMIS:ZeroPosterior','Posterior resampling requires positive finite total mass.');
logMass = logMass - normalizer;
end
