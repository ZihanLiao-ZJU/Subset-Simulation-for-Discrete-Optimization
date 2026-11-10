function [CoV_lw,CoV_up,parts] = MisCov(logEstimate,logW,logAlpha,logH,chainLengths,normalizerJacobian,independentSamples,sampleDesign)
% Reliability-v1 CoV construction for any nonnegative target f(u).
% Reliability uses f=I(g<0)*p; evidence uses f=L*p. logW{i} contains the
% balance-weighted IS contributions whose level means sum to the estimate.
% Retains v1's MCMC, normalizer and cross terms (Eqns. 34-58, Appendix B/C).
% Terms are divided by Z^2. Derivatives w.r.t. log C_j cancel C_j exactly,
% avoiding exponentiation of tiny normalizers without changing the formulas.
Nite = numel(logW);
if nargin<8||isempty(sampleDesign),sampleDesign=cell(Nite,1);end
assert(iscell(sampleDesign)&&numel(sampleDesign)==Nite,'SeMIS:SampleDesign','One design per sampled level is required.');
if nargin<5 || isempty(chainLengths), chainLengths=cell(Nite,1); end
assert(iscell(chainLengths) && numel(chainLengths)==Nite, ...
    'SeMIS:ChainLayout','Provide one optional chain-length vector per level.');
if nargin<7 || isempty(independentSamples)
    independentSamples=false(Nite,1);
    if Nite>0, independentSamples(1)=true; end
end
independentSamples=MISIndependentSamples(struct('IndependentSamples',independentSamples, ...
    'ChainLengths',{chainLengths}),Nite);
M = Nite-1;
if nargin<6, normalizerJacobian=[]; end
if ~isempty(normalizerJacobian)
    assert(isnumeric(normalizerJacobian) && isreal(normalizerJacobian) && ...
        isequal(size(normalizerJacobian),[Nite,M]) && ...
        all(isfinite(normalizerJacobian),'all') && ...
        all(normalizerJacobian==0 | normalizerJacobian==1,'all') && ...
        all(normalizerJacobian(1,:)==0), ...
        'SeMIS:NormalizerJacobian','Provide a binary [levels,transitions] normalizer Jacobian.');
    assert(all(triu(normalizerJacobian(2:end,:),1)==0,'all'), ...
        'SeMIS:NormalizerJacobian','Normalizers cannot depend on later transitions.');
end
assert(numel(logAlpha)==Nite && numel(logH)==M, ...
    'SeMIS:MISSize','CoV histories have inconsistent lengths.');
if any(~cellfun(@isempty,sampleDesign))
    assert(M==0||(~isempty(normalizerJacobian)&&~any(normalizerJacobian,'all')), ...
        'SeMIS:StratifiedNormalizer','Stratified variance currently requires known absolute normalizers.');
end
parts = struct('mcmc',zeros(Nite,1),'normalizerLower',zeros(Nite,1), ...
    'normalizerUpper',zeros(Nite,1),'crossUpper',zeros(Nite,1));
if logEstimate==-inf
    CoV_lw = inf; CoV_up = inf; return
end
assert(isreal(logEstimate) && isfinite(logEstimate), ...
    'SeMIS:InvalidEstimate','The estimate must have a real finite logarithm.');
delta2 = zeros(M,1);
for m = 1:M
    if ~isempty(normalizerJacobian) && ~any(normalizerJacobian(:,m)), continue; end
    h = logH{m}-logmean(logH{m},'all');
    if independentSamples(m)
        delta2(m) = var(exp(h(:)),0)/numel(h);
    else
        delta2(m) = MoveBlock(h,chainLengths{m});
    end
end
delta = sqrt(max(delta2,0));
s = cumsum(delta);
s2 = cumsum(delta2);
CovLower = zeros(M,M);
for j = 1:M
    for l = 1:M
        CovLower(j,l) = s2(min(j,l));
    end
end
CovUpper = s*s.';
biasRelative = max(0.5*(s.^2-s2),0);
if ~isempty(normalizerJacobian)
    B=normalizerJacobian(2:end,:);
    CovLower=(B.*delta2.')*B.';
    cumulative=B*delta;
    CovUpper=cumulative*cumulative.';
    biasRelative=max(0.5*(cumulative.^2-B*delta2),0);
end
for i = 1:Nite
    normalizedLogW = logW{i}-logEstimate;
    W = exp(normalizedLogW(:));
    if ~isempty(sampleDesign{i})
        design=MISSampleDesign(sampleDesign{i},numel(W));v=0;
        for k=1:numel(design.counts)
            v=v+design.counts(k)*var(W(design.labels==k),0);
        end
        parts.mcmc(i)=v/numel(W)^2;
    elseif independentSamples(i)
        parts.mcmc(i) = var(W,0)/numel(W);
    else
        parts.mcmc(i) = MoveBlock(normalizedLogW,chainLengths{i});
    end
    if M==0, continue; end
    alpha = exp(logAlpha{i}(2:end,:)).';
    % A_j = C_j/Z * d P_i/d C_j = mean((W/Z)*alpha_j).
    dW = W.*alpha;
    A = mean(dW,1).';
    parts.normalizerLower(i) = max(A.'*CovLower*A,0);
    parts.normalizerUpper(i) = max(A.'*CovUpper*A,0);
    centered = W-mean(W);
    ss = sum(centered.^2);
    if ss>0
        scale = parts.mcmc(i)/ss;
        F = (2*scale*sum(centered.*(dW-mean(dW,1)),1)).';
        % The bound concerns |bias(C_j)|, not its sign: |F| is required
        % for an upper bound rather than allowing a negative cross correction.
        parts.crossUpper(i) = sum(abs(F).*biasRelative);
    end
end
CoV_lw = sqrt(sum(parts.mcmc+parts.normalizerLower));
CoV_up = sum(sqrt(parts.mcmc+parts.normalizerUpper+parts.crossUpper));
% These are v1 plug-in correlation bounds under its assumptions, including
% nonnegative between-level covariance for the lower estimate.
end
