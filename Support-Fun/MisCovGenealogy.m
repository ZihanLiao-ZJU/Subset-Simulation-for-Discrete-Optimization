function report=MisCovGenealogy(context,logEstimate,logW,logAlpha,logH,logC,B,options) %#ok<INUSD>
% Joint MIS influence scores clustered by the actual initial ancestor.
% Approximate root independence; NOT an unbiased SMC/splitting theorem.
% All initial roots remain, even if extinct. No calls, resampling or blocks.
if nargin<8,options=struct();end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),{'FixedDesign'})), ...
    'SeMIS:UncertaintyOptions','Only FixedDesign is supported by the genealogy backend.');
fixed=false;
if isfield(options,'FixedDesign')
    assert(islogical(options.FixedDesign)&&isscalar(options.FixedDesign),'SeMIS:UncertaintyOptions','FixedDesign must be logical.');
    fixed=options.FixedDesign;
end
report=struct('method','joint-first-order-root-genealogy','scope','fixed-realised-design delta and root-cluster approximation', ...
    'status','approximate','cov',NaN,'relativeVariance',NaN,'variance',NaN,'logVariance',NaN, ...
    'logEstimate',logEstimate,'fixedDesignDeclared',fixed,'fullAdaptiveVarianceValidated',false, ...
    'limitations',{{'Initial-root clusters are treated as approximately independent; resampling can couple distinct roots even under fixed ISDs.', ...
    'The initial-root Bessel factor is a cluster correction, not a general resampling-unbiasedness result.', ...
    'Lineage depletion can cause loss of late-level variance, not just noisy estimates.'}});
if isfield(context,'SeedInfo') && any(cellfun(@(s)isstruct(s) && ...
        isfield(s,'selection') && strcmp(s.selection,'custom'),context.SeedInfo))
    report.status='unsupported_seed_selection';
    report.limitations{end+1}=['The recorded custom parent-selection law has no implemented selection-variance ', ...
        'model in this backend. Point MIS integration remains available; covariance is not identified here.'];
    return
end
if ~isfinite(logEstimate),report.status='zero_or_nonfinite_estimate';return;end
L=numel(logW);N=cellfun(@numel,logW);N0=N(1);
if N0<2||~context.IndependentSamples(1)||~isempty(context.SampleDesign{1})
    report.status='unsupported_initial_root_design';return
end
if context.Initial.cov~=0
    report.status='unsupported_external_initial_uncertainty';return
end
% The first implementation requires one connected ancestral population.
% Unrelated fresh populations and fitted direct proposals need another model.
for i=2:L
    if context.IndependentSamples(i)||~isempty(context.SampleDesign{i})
        report.status='unsupported_unconnected_population';return
    end
    if ~isfield(context,'SeedInfo')||numel(context.SeedInfo)<i||isempty(context.SeedInfo{i})
        report.status='missing_seed_ancestry';return
    end
    seed=context.SeedInfo{i};
    if ~isfield(seed,'integrationPreserved')||~isequal(seed.integrationPreserved,true)
        report.status='unverified_seed_lift_transport';return
    end
end
[W,R,psi,A,gamma,~,status]=MISJointInfluence(context,logEstimate,logW,logAlpha,logH,B);
if ~isempty(status),report.status=status;return;end
factor=N0/(N0-1);RS=zeros(N0,L);RN=RS;rootCounts=RS;
same=zeros(L,3);chains=same;counts=zeros(L,1);roots=(1:N0).';
for i=1:L
    if i==1
        groups=num2cell((1:N0).');
    else
        groups=chainGroups(context,i,N(i));p=context.SeedInfo{i}.parentIndices(:);
        assert(numel(p)==numel(groups)&&all(isfinite(p)&p==floor(p)&p>=1&p<=N(i-1)), ...
            'SeMIS:SeedHistory','Actual parent indices must identify one parent per retained chain.');
        parentRoots=roots(p);roots=zeros(N(i),1);
        for c=1:numel(groups),roots(groups{c})=parentRoots(c);end
    end
    assert(all(roots>=1&roots<=N0),'SeMIS:SeedHistory','Incomplete root assignment.');
    v=[W{i},gamma(i)*R{i}];v=(v-mean(v,1))/N(i);
    if any(~isfinite(v),'all'),report.status='nonfinite_influence';return;end
    RS(:,i)=accumarray(roots,v(:,1),[N0,1]);RN(:,i)=accumarray(roots,v(:,2),[N0,1]);
    rootCounts(:,i)=accumarray(roots,ones(N(i),1),[N0,1]);
    same(i,:)=factor*[sum(v(:,1).^2),sum(v(:,2).^2),2*sum(v(:,1).*v(:,2))];
    cs=zeros(numel(groups),2);
    for c=1:numel(groups),cs(c,:)=sum(v(groups{c},:),1);end
    chains(i,:)=factor*[sum(cs(:,1).^2),sum(cs(:,2).^2),2*sum(cs(:,1).*cs(:,2))];
    counts(i)=numel(groups);
    % The joint score agrees with separately centred source observables.
    joint=accumarray(roots,(psi{i}-mean(psi{i}))/N(i),[N0,1]);
    assert(max(abs(joint-RS(:,i)-RN(:,i)))<1e-10*max(1,max(abs(joint))));
end
scores=RS+RN;C=factor*(scores.'*scores);total=factor*sum(sum(scores,2).^2);
assert(abs(sum(C,'all')-total)<1e-10*max(1,total));
levelSource=factor*[sum(RS.^2,1).',sum(RN.^2,1).',2*sum(RS.*RN,1).'];
allS=sum(RS,2);allN=sum(RN,2);
sourceTotal=factor*[sum(allS.^2),sum(allN.^2),2*sum(allS.*allN)];
partition=[sum(same,1);sum(chains-same,1);sum(levelSource-chains,1);sourceTotal-sum(levelSource,1)];
assert(abs(sum(partition,'all')-total)<1e-10*max(1,total));
report.relativeVariance=total;report.cov=sqrt(total);
if total==0,report.logVariance=-inf;report.variance=0;
else,report.logVariance=2*logEstimate+log(total);report.variance=exp(report.logVariance);end
report.covariance=C;report.withinLevel=diag(C);report.rootScores=scores;
report.rootSamplingScores=RS;report.rootNormalizerScores=RN;report.rootCounts=rootCounts;
report.rootCount=N0;report.rootAlive=sum(rootCounts>0,1).';
report.rootEffective=(N(:).^2)./sum(rootCounts.^2,1).';
report.rootLargestFraction=max(rootCounts,[],1).'./N(:);
report.rootCollapsed=report.rootAlive<=1;report.rootCorrection=factor;report.chainCount=counts;
report.sensitivity=A;report.incrementSensitivity=gamma(1:end-1);
report.sourceDependence=partition;
report.sourceLabels={'sampling','normalizer','signedCross'};
report.dependenceLabels={'sameState','intrachain','interchain','interlevel'};
report.components=struct('withinSampling',sum(levelSource(:,1)), ...
    'withinNormalizer',sum(levelSource(:,2)),'withinSamplingNormalizer',sum(levelSource(:,3)), ...
    'interlevel',sum(partition(4,:)),'initial',0);
if any(report.rootCollapsed(2:end))
    report.status='lineage_degenerate';
end
if ~fixed
    report.limitations{end+1}='Adaptive thresholds, fitted proposals, shared sampler adaptation and random stopping are held at observed values; their unconditional effects are not fully identified.';
end
end

function groups=chainGroups(context,i,N)
lengths=context.ChainLengths{i};
if ~isempty(lengths)
    assert(all(lengths>0&lengths==floor(lengths))&&sum(lengths)==N,'SeMIS:ChainLayout','Invalid chain lengths.');
    edge=[0;cumsum(lengths(:))];groups=cell(numel(lengths),1);
    for c=1:numel(lengths),groups{c}=(edge(c)+1:edge(c+1)).';end
else
    shape=size(context.Yraw{i},2:3);nc=shape(1);ns=shape(2);
    assert(nc*ns==N,'SeMIS:ChainLayout','Invalid rectangular chain layout.');
    groups=cell(nc,1);for c=1:nc,groups{c}=(c:nc:nc*ns).';end
end
end
