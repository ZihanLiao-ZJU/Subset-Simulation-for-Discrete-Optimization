function independent = MISIndependentSamples(in,Nite)
% Per-level conditional iid flag. Legacy histories assume iid only initially.
independent=false(Nite,1);
if Nite>0, independent(1)=true; end
if ~isfield(in,'IndependentSamples'), return; end
value=in.IndependentSamples;
if ~islogical(value) || ~isvector(value) || numel(value)~=Nite
    error('SeMIS:IndependentSamples','Provide one logical independent-sampling flag per level.');
end
independent=value(:);
if isfield(in,'ChainLengths')
    if ~iscell(in.ChainLengths) || numel(in.ChainLengths)~=Nite
        error('SeMIS:ChainLayout','Provide one optional chain-length vector per level.');
    end
    for level=find(independent).'
        if ~isempty(in.ChainLengths{level}) && any(in.ChainLengths{level}~=1)
            error('SeMIS:IndependentSamples','Independent samples cannot carry multi-step chain lengths.');
        end
    end
end
end
