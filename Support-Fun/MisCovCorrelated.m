function report=MisCovCorrelated(context,logEstimate,logW,logAlpha,logH,logC,B,options)
% Joint first-order MIS uncertainty with Miao--Low inspired dependence.
% The paper establishes a model for product-SuS. Continuous weighted seed
% transport, ragged chains and MIS influence functions below are extensions,
% NOT a proof for arbitrary adaptive proposals or stopping rules.
% No model evaluations or random draws are made by this postprocessor.
if nargin<8, options=struct(); end
assert(isstruct(options)&&isscalar(options),'SeMIS:UncertaintyOptions','Options must be a scalar struct.');
allowed={'FixedDesign'};
assert(all(ismember(fieldnames(options),allowed)), ...
    'SeMIS:UncertaintyOptions','Unknown uncertainty option.');
fixed=false;
if isfield(options,'FixedDesign')
    assert(islogical(options.FixedDesign)&&isscalar(options.FixedDesign), ...
        'SeMIS:UncertaintyOptions','FixedDesign must be logical.');
    fixed=options.FixedDesign;
end
L=numel(logW); N=cellfun(@numel,logW);
report=struct('method','joint-first-order-block-seed', ...
    'reference','Miao and Low (2026), doi:10.1016/j.strusafe.2026.102690', ...
    'scope','fixed-realised-design delta approximation with seed mediation', ...
    'status','approximate','cov',NaN,'relativeVariance',NaN,'variance',NaN, ...
    'logVariance',NaN,'logEstimate',logEstimate,'limitations',{{}}, ...
    'fullAdaptiveVarianceValidated',false,'fixedDesignDeclared',fixed);
if isfield(context,'SeedInfo') && any(cellfun(@(s)isstruct(s) && ...
        isfield(s,'selection') && strcmp(s.selection,'custom'),context.SeedInfo))
    report.status='unsupported_seed_selection';
    report.limitations={['The recorded custom parent-selection law has no implemented seed-selection ', ...
        'variance model. Point MIS integration remains available; covariance is not identified here.']};
    return
end
if ~isfinite(logEstimate)
    report.status='zero_or_nonfinite_estimate';
    report.limitations={'A zero/nonfinite estimate does not identify its sampling variance.'};
    return
end
[W,R,psi,A,gamma,logHmean,influenceStatus]=MISJointInfluence(context,logEstimate,logW,logAlpha,logH,B);
if ~isempty(influenceStatus)
    report.status=influenceStatus;
    report.limitations={'A required transition ratio has no finite positive empirical overlap.'};
    return
end
% A canonical ancestry order is reconstructed for analysis only. The actual
% seed permutation, sampler adaptation order and simulation are unchanged.
groups=cell(L,1); order=cell(L,1); ranks=cell(L,1);
within=zeros(L,1); withinSampling=zeros(L,1); withinNormalizer=zeros(L,1);
withinCross=zeros(L,1); blocks=zeros(L,1); levelInfo=cell(L,1);
for i=1:L
    groups{i}=chainGroups(context,i,N(i));
    nc=numel(groups{i}); order{i}=1:nc;
    seed=seedAt(context,i);
    if i>1 && ~isempty(seed)
        parent=seed.parentIndices(:);
        assert(numel(parent)==nc,'SeMIS:SeedHistory','One parent is required per retained chain.');
        [~,order{i}]=sortrows([ranks{i-1}(parent), (1:nc).'],[1,2]);
        order{i}=order{i}.';
    elseif i>1 && ~context.IndependentSamples(i) && isempty(context.SampleDesign{i})
        report.limitations{end+1}=sprintf('Level %d has no recorded seed ancestry; ordered interchain dependence is not identified.',i);
    end
    ranks{i}=zeros(N(i),1); rank=0;
    for c=order{i}
        ix=groups{i}{c}; ranks{i}(ix)=rank+(1:numel(ix)); rank=rank+numel(ix);
    end
    obs=[W{i},gamma(i)*R{i},psi{i}];
    [S,b,detail]=withinStats(obs,groups{i},order{i},context,i);
    within(i)=S(3,3); withinSampling(i)=S(1,1);
    withinNormalizer(i)=S(2,2); withinCross(i)=2*S(1,2);
    blocks(i)=b; levelInfo{i}=detail;
    if nc==1 && N(i)>1 && ~context.IndependentSamples(i)
        report.limitations{end+1}=sprintf('Level %d has one retained chain; temporal blocks replace across-chain composites.',i);
    end
    if numel(unique(cellfun(@numel,groups{i})))>1
        report.limitations{end+1}=sprintf('Level %d uses length-adjusted ragged-chain composites, an extension beyond equal-chain SuS.',i);
    end
end
if any(~isfinite(within))
    report.status='insufficient_variance_information';
    report.limitations{end+1}='At least one level has fewer than two effective observations or nonfinite moments.';
    return
end

adj=zeros(max(0,L-1),1); adjacent=cell(max(0,L-1),1);
for i=1:L-1
    j=i+1; seed=seedAt(context,j);
    d=struct('sourceLevel',i,'childLevel',j,'rho',0,'levelToSeed',0, ...
        'seedToChain',0,'selectionVariance',0,'status','');
    if isempty(seed)
        if context.IndependentSamples(j)||~isempty(context.SampleDesign{j})
            d.status='direct_no_seed_mediation';
            if ~fixed
                report.limitations{end+1}=sprintf('Transition %d->%d: direct adaptive-proposal learning covariance is not estimated.',i,j);
            end
        else
            d.status='missing_seed_ancestry';
        end
        adjacent{i}=d; continue
    end
    if ~isfield(seed,'integrationPreserved') || ~isequal(seed.integrationPreserved,true)
        d.status='unverified_seed_lift_transport';
        report.limitations{end+1}=sprintf('Transition %d->%d has no certified lift-invariant integration observable.',i,j);
        adjacent{i}=d; continue
    end
    % Evaluate the child influence on parent records using cached balance
    % densities. Zero transition support is masked before any quotient.
    active=isfinite(logAlpha{i}(j,:)).';
    r=R{i}; f=zeros(N(i),1);
    f(active)=N(j)/N(i)*W{i}(active);
    if gamma(j)~=0
        logRatio=logAlpha{i}(j+1,active)-logAlpha{i}(j,active) ...
            +log(N(j)/N(j+1))+logC(j+1)-logC(j)-logHmean(j);
        f(active)=f(active)+gamma(j)*exp(logRatio(:));
    end
    if any(~isfinite(f(active))) || any(~isfinite(r))
        d.status='nonfinite_seed_transport';
        report.limitations{end+1}=sprintf('Transition %d->%d seed transport has nonfinite moments.',i,j);
        adjacent{i}=d; continue
    end
    % Density support and actual selection support can differ at boundary
    % ties. Normalizer influences retain R; seed transport uses the design
    % that was actually used to select parents.
    transport=r;
    if strcmp(seed.selection,'survivors')
        if ~isfield(seed,'candidateIndices')
            d.status='missing_survivor_pool';
            adjacent{i}=d;
            report.limitations{end+1}=sprintf('Transition %d->%d lacks the actual eligible seed pool.',i,j);
            continue
        end
        pool=seed.candidateIndices(:);
        transport=zeros(size(r)); transport(pool)=N(i)/numel(pool);
    end
    p=transport/sum(transport); T=sum(p.*f); v=transport.*(f-T);
    parent=seed.parentIndices(:);
    if any(~active(parent))
        error('SeMIS:SeedHistory','Recorded parent lies outside child proposal support.');
    end
    lengths=cellfun(@numel,groups{j}); omega=lengths(:)/sum(lengths);
    if strcmp(seed.selection,'weighted')
        selectionVar=sum(omega.^2)*sum(p.*(f-T).^2);
    elseif strcmp(seed.selection,'survivors')
        % The actual seed pool can exclude boundary ties. Preserve its
        % indices; do not infer that pool from density support after the fact.
        K=numel(pool);
        if K<2, selectionVar=0;
        else, selectionVar=max(0,sum(omega.^2)-1/K)*var(f(pool),0); end
    else
        d.status='unsupported_seed_selection';
        report.limitations{end+1}=sprintf('Transition %d->%d has an unsupported seed-selection design.',i,j);
        adjacent{i}=d; continue
    end
    % The ratio's variability excludes randomness in the realised resampling
    % step. Add that design variance before correlating with actual seeds.
    C=composites([psi{i},v],groups{i});
    C=C(order{i},:);
    mode='ordered-blocks';
    if context.IndependentSamples(i)||firstGenerationIndependent(context,i)
        mode='independent';
    end
    pair=MiaoLowBlockStats(C,struct('mode',mode));
    if ~pair.available
        d.status='unidentifiable_parent_composites'; adjacent{i}=d; continue
    end
    first=nan(numel(pair.jointBValues),1);
    for b=1:numel(first)
        V=pair.jointCovByBlock(:,:,b);
        denominator=sqrt(max(V(1,1),0)*max(V(2,2)+selectionVar,0));
        if denominator>0, first(b)=V(1,2)/denominator; end
    end
    % Max over b is the paper's positive-dependence selection rule. Negative
    % values are retained if all candidates are negative, not turned into |rho|.
    [rho1,k1]=max(first,[],'omitnan');
    sourceSeed=f(parent);
    mu=sum(omega.*sourceSeed);
    seedComposite=numel(omega)*omega.*(sourceSeed-mu);
    childComposite=composites(psi{j},groups{j});
    childPair=[seedComposite,childComposite]; childPair=childPair(order{j},:);
    mode='ordered-blocks';
    if firstGenerationIndependent(context,j), mode='independent'; end
    next=MiaoLowBlockStats(childPair,struct('mode',mode));
    if ~next.available
        d.status='unidentifiable_child_composites'; adjacent{i}=d; continue
    end
    second=nan(numel(next.jointBValues),1);
    for b=1:numel(second)
        V=next.jointCovByBlock(:,:,b);
        denominator=sqrt(max(V(1,1),0)*max(V(2,2),0));
        if denominator>0, second(b)=V(1,2)/denominator; end
    end
    [rho2,k2]=max(second,[],'omitnan');
    if ~isfinite(rho1)||~isfinite(rho2)
        d.status='zero_empirical_seed_variance';
        report.limitations{end+1}=sprintf('Transition %d->%d has a zero empirical variance in a seed-correlation leg; zero covariance is a placeholder, not identified independence.',i,j);
        adjacent{i}=d; continue
    end
    assert(abs(rho1)<=1+1e-10 && abs(rho2)<=1+1e-10, ...
        'SeMIS:CorrelationRange','Estimated seed correlation is outside [-1,1].');
    rho1=max(-1,min(1,rho1)); rho2=max(-1,min(1,rho2));
    adj(i)=rho1*rho2;
    d.rho=adj(i); d.levelToSeed=rho1; d.seedToChain=rho2;
    d.selectionVariance=selectionVar; d.levelSeedBlock=pair.jointBValues(k1);
    d.seedChainBlock=next.jointBValues(k2); d.status='seed_mediation_approximation';
    adjacent{i}=d;
end
correlation=eye(L);
for i=1:L
    for j=i+1:L
        correlation(i,j)=prod(adj(i:j-1)); correlation(j,i)=correlation(i,j);
    end
end
covariance=correlation.*(sqrt(max(within,0))*sqrt(max(within,0)).');
conditional=sum(covariance,'all');
assert(conditional>=-1e-10*max(1,sum(within)), ...
    'SeMIS:NegativeVariance','The correlation model produced negative variance.');
conditional=max(0,conditional);
initial=context.Initial.cov^2;
total=initial+(1+initial)*conditional;
report.relativeVariance=total; report.cov=sqrt(total);
if total==0
    report.logVariance=-inf; report.variance=0;
else
    report.logVariance=2*logEstimate+log(total);
    report.variance=exp(report.logVariance);
end
report.covariance=covariance; report.correlation=correlation;
report.sensitivity=A; report.incrementSensitivity=gamma(1:end-1);
report.withinLevel=within; report.withinBlock=blocks; report.levelInfo=levelInfo;
report.adjacent=adjacent;
report.components=struct('withinSampling',sum(withinSampling), ...
    'withinNormalizer',sum(withinNormalizer), ...
    'withinSamplingNormalizer',sum(withinCross), ...
    'interlevel',sum(covariance,'all')-sum(within), ...
    'initial',initial*(1+conditional));
report.conditionalRelativeVariance=conditional;
if ~fixed
    report.limitations{end+1}=['Proposal fitting, adaptive thresholds, shared sampler-scale adaptation and random stopping ' ...
        'are held at their realised values; the single-run approximation is not a general unconditional variance guarantee.'];
end
if L>1
    report.limitations{end+1}=['Adjacent scalar seed mediation and products of adjacent correlations ' ...
        'are model closures, not consequences of the MIS delta expansion.'];
end
incomplete=cellfun(@(d)~strcmp(d.status,'seed_mediation_approximation') && ...
    ~(fixed&&strcmp(d.status,'direct_no_seed_mediation')),adjacent);
if any(incomplete)
    report.status='partial_dependence_model';
end
report.limitations=unique(report.limitations,'stable');
end

function seed=seedAt(context,i)
seed=[];
if isfield(context,'SeedInfo') && numel(context.SeedInfo)>=i
    seed=context.SeedInfo{i};
end
end

function yes=firstGenerationIndependent(context,i)
% Paper convention conditional on the realised threshold and transition
% kernel: distinct survivors from an iid parent initiate independent chains.
% Weighted resampling (even with incidentally distinct indices) does not.
yes=false;
if i<=1 || ~context.IndependentSamples(i-1), return; end
seed=seedAt(context,i);
if isempty(seed)||~isfield(seed,'selection')|| ...
        ~strcmp(seed.selection,'survivors')|| ...
        ~isfield(seed,'integrationPreserved')||~isequal(seed.integrationPreserved,true)
    return
end
parents=seed.parentIndices(:);
yes=~isempty(parents)&&numel(unique(parents))==numel(parents);
end

function groups=chainGroups(context,i,N)
lengths=context.ChainLengths{i};
if context.IndependentSamples(i)||~isempty(context.SampleDesign{i})
    groups=num2cell((1:N).'); return
end
if ~isempty(lengths)
    edge=[0;cumsum(lengths(:))]; groups=cell(numel(lengths),1);
    for c=1:numel(lengths), groups{c}=(edge(c)+1:edge(c+1)).'; end
else
    shape=size(context.Yraw{i},2:3); nc=shape(1); ns=shape(2);
    groups=cell(nc,1);
    for c=1:nc, groups{c}=(c:nc:nc*ns).'; end
end
end

function C=composites(obs,groups)
N=size(obs,1); nc=numel(groups); mu=mean(obs,1);
C=zeros(nc,size(obs,2));
for c=1:nc
    ix=groups{c};
    % Length-adjusted centered totals recover chain means up to a constant
    % for equal lengths and are zero for a constant observable at all lengths.
    C(c,:)=nc/N*sum(obs(ix,:)-mu,1);
end
end

function [S,b,detail]=withinStats(obs,groups,order,context,i)
N=size(obs,1); P=size(obs,2); design=context.SampleDesign{i};
detail=struct('mode','','chainCount',numel(groups));
if ~isempty(design)
    S=zeros(P); b=1; detail.mode='fixed-stratum';
    for k=1:numel(design.counts)
        if design.counts(k)<2
            S=nan(P); detail.mode='insufficient_stratum_observations'; return
        end
        S=S+design.counts(k)*cov(obs(design.labels==k,:),0)/N^2;
    end
    return
end
if context.IndependentSamples(i)
    if N<2, S=nan(P); else, S=cov(obs,0)/N; end
    b=1; detail.mode='conditional-iid'; return
end
if numel(groups)==1
    C=obs-mean(obs,1); detail.mode='single-chain-temporal-blocks';
else
    C=composites(obs,groups); C=C(order,:); detail.mode='ancestry-ordered-chain-blocks';
end
mode='ordered-blocks';
if numel(groups)>1 && firstGenerationIndependent(context,i)
    mode='independent'; detail.mode='first-generation-independent-chains';
end
s=MiaoLowBlockStats(C,struct('mode',mode));
[~,k]=max(squeeze(s.covByBlock(3,3,:)));
S=s.covByBlock(:,:,k); b=s.bValues(k);
end
