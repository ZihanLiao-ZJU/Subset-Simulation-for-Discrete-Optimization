function report=MISUncertainty(context,logEstimate,logW,logAlpha,logH,logC,B,options)
% Select the covariance estimator without changing MIS or its influence map.
if nargin<8,options=struct();end
assert(isstruct(options)&&isscalar(options),'SeMIS:UncertaintyOptions','Options must be a scalar struct.');
backend='miao';
if isfield(options,'Backend')
    backend=options.Backend;options=rmfield(options,'Backend');
end
assert((ischar(backend)||isstring(backend))&&isscalar(string(backend))&& ...
    any(strcmp(backend,{'miao','genealogy'})),'SeMIS:UncertaintyBackend','Backend must be miao or genealogy.');
if strcmp(backend,'miao')
    report=MisCovCorrelated(context,logEstimate,logW,logAlpha,logH,logC,B,options);
else
    report=MisCovGenealogy(context,logEstimate,logW,logAlpha,logH,logC,B,options);
end
end
