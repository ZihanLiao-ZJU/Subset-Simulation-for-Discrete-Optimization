function stats = MiaoLowBlockStats(composites,options)
%MIAOLOWBLOCKSTATS Joint covariance of ordered, equal-weight chain composites.
% Each row is one chain composite, each column a REAL SIGNED observable.
% The caller must construct composites in documented ancestral order and
% handle unequal chain lengths or fixed-allocation weights before calling.
%
% Miao & Low (2026), Structural Safety 120, 102690, Eqs. (9)-(16),
% (50)-(53), (64)-(66): marginal blocks are capped at floor(Nc/8);
% correlation candidates at floor(Nc/20). No block length is optimized
% separately for covariance entries: each returned joint matrix is PSD
% up to floating-point error. Pair/ratio correlation selection is external.
%
% Explicit b=1 convention: unbiased covariance(composites)/Nc. For b>=2,
% covariance(block means, population normalization)*b/Nc follows Eq. (14).
% Including b=1 as a candidate and falling back to it below the paper's
% minimum chain counts are implementation conventions, recorded below.
% This helper does not establish independence, infer ancestry, estimate
% interlevel correlations, or change the legacy MoveBlock implementation.

if nargin<2 || isempty(options), options=struct; end
assert(isnumeric(composites)&&isreal(composites)&&ismatrix(composites) && ...
    size(composites,2)>=1&&all(isfinite(composites),'all'), ...
    'SeMIS:MiaoLowObservables','Provide a finite real Nc-by-P matrix.');
assert(isstruct(options)&&isscalar(options), ...
    'SeMIS:MiaoLowOptions','Options must be a scalar structure.');
assert(all(ismember(fieldnames(options),{'mode'})), ...
    'SeMIS:MiaoLowOptions','The only supported option is mode.');
if ~isfield(options,'mode'), options.mode='ordered-blocks'; end
mode=validatestring(options.mode,{'independent','ordered-blocks'});
C=double(composites);
[nc,p]=size(C);
capVar=floor(nc/8); capJoint=floor(nc/20);
if strcmp(mode,'independent')
    bValues=1;
else
    bValues=1:max(1,capVar);
end
jointMask=bValues<=max(1,capJoint);
nb=numel(bValues);
stats=struct('available',false,'status','insufficient_chains', ...
    'mode',mode,'nChains',nc,'nObservables',p, ...
    'compositeCovariance',nan(p),'baseCovariance',nan(p), ...
    'bValues',bValues,'covByBlock',nan(p,p,nb), ...
    'correlationByBlock',nan(p,p,nb), ...
    'varianceMax',nan(1,p),'selectedBlock',nan(1,p), ...
    'jointMask',jointMask,'jointBValues',bValues(jointMask), ...
    'jointCovByBlock',nan(p,p,nnz(jointMask)), ...
    'jointCorrelationByBlock',nan(p,p,nnz(jointMask)), ...
    'paperVarianceBlockLimit',capVar,'paperCorrelationBlockLimit',capJoint, ...
    'marginalBlockFallback',strcmp(mode,'ordered-blocks')&&capVar<1, ...
    'correlationBlockFallback',strcmp(mode,'ordered-blocks')&&capJoint<1, ...
    'b1Convention','Unbiased chain-composite covariance divided by Nc; included as baseline.', ...
    'orderingAssumed',strcmp(mode,'ordered-blocks'));
if nc<2, return; end

centered=C-mean(C,1);
chainCov=(centered.'*centered)/(nc-1);
if any(~isfinite(chainCov),'all')
    stats.status='numerical_overflow'; return
end
stats.compositeCovariance=(chainCov+chainCov.')/2;
stats.baseCovariance=stats.compositeCovariance/nc;
covariances=nan(p,p,nb);
covariances(:,:,1)=stats.baseCovariance;
% Removing the global mean before cumulative summation improves stability.
sums=[zeros(1,p);cumsum(centered,1)];
for k=2:nb
    b=bValues(k); nm=nc-b+1;
    blocks=(sums(b+1:end,:)-sums(1:end-b,:))/b;
    blocks=blocks-mean(blocks,1);
    V=(blocks.'*blocks)*(b/(nm*nc));
    if any(~isfinite(V),'all')
        stats.status='numerical_overflow'; return
    end
    covariances(:,:,k)=(V+V.')/2;
end
correlations=nan(p,p,nb);
diagonals=zeros(nb,p);
for k=1:nb
    V=covariances(:,:,k); v=diag(V);
    diagonals(k,:)=v.';
    positive=find(v>0);
    if ~isempty(positive)
        sd=sqrt(v(positive));
        R=(V(positive,positive)./sd)./sd.';
        % A material violation is an error, not a truncation of evidence.
        tolerance=64*eps;
        assert(all(abs(R)<=1+tolerance,'all'), ...
            'SeMIS:MiaoLowCorrelation','Joint covariance gives an invalid correlation.');
        R(abs(R)>1)=sign(R(abs(R)>1));
        correlations(positive,positive,k)=R;
    end
end
[stats.varianceMax,index]=max(diagonals,[],1);
stats.selectedBlock=bValues(index);
stats.covByBlock=covariances;
stats.correlationByBlock=correlations;
stats.jointCovByBlock=covariances(:,:,jointMask);
stats.jointCorrelationByBlock=correlations(:,:,jointMask);
stats.available=true;
stats.status='ok';
end
