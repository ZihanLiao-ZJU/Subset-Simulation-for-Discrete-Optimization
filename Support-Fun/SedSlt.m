function [x_sed,y_sed,intBay,Ncal,seedInfo] = SedSlt(intBay,x,y,Nsam,fixedChains,selection)
% Resample in integration space, then lift seeds into sampling space.
% Identity ISDs retain the original rectangular-chain resampling behavior.
validateattributes(Nsam,{'numeric'},{'scalar','integer','positive','finite'});
fixed = nargin>=5 && ~isempty(fixedChains);
if nargin<6, selection='auto'; end
selection=validatestring(selection,{'auto','weighted','survivors'});
if fixed, validateattributes(fixedChains,{'numeric'},{'scalar','integer','positive','finite'}); end
Ncal = 0;
seedInfo = [];
x = reshape(x,size(x,1),[]);
y = reshape(y,size(y,1),[]);
view = ISDIntegrationView(intBay);
logw = ResRat(view,x,y,intBay.G);
w = exp(logw(:));
if ~any(w>0)
    x_sed = zeros(size(x,1),0); y_sed = zeros(size(y,1),0); return
end
% Retained mass controls chain count. A divisor of the requested budget
% avoids silently changing Nsam by rounding the number of steps per chain.
if fixed
    Nc = fixedChains;
else
    mass = sum(w);
    % A solved retention constraint can land a few ulps below an integer.
    % Flooring 99.99999999999997 would halve a 1000-point layout from
    % 100x10 to 50x20. Snap roundoff only; genuine fractional mass still floors.
    nearest = round(mass);
    if abs(mass-nearest)<=32*eps(max(1,mass)), mass=nearest; end
    desired = max(1,min(Nsam,floor(mass)));
    d = 1:floor(sqrt(Nsam));
    d = d(mod(Nsam,d)==0);
    divisors = unique([d,Nsam./d]);
    Nc = max(divisors(divisors<=desired));
end
if strcmp(selection,'survivors')
    % Re-evaluating a Gaussian in a different batch layout can change its
    % log density by a few ulps. A hard event still has constant finite
    % transition weights; use a cancellation-scale tolerance for that test.
    densityScale=y(1:2,:); densityScale=densityScale(isfinite(densityScale));
    tolerance=32*eps(max([1;abs(densityScale(:))]));
    active=isfinite(logw(:));
    binary=all(abs(logw(active))<=tolerance);
    if ismethod(intBay,'SeedEligibility')
        % An ISD may exclude finite-history boundary ties from its seed
        % pool without changing the density used by MIS or MH proposals.
        eligible=intBay.SeedEligibility(x,y);
        if ~islogical(eligible) || ~isvector(eligible) || ...
                numel(eligible)~=numel(active) || any(eligible(:) & ~active)
            error('SeMIS:SeedEligibility', ...
                'SeedEligibility must return one logical per point, selecting only current ISD support.');
        end
        active=active & eligible(:);
    end
    if ~binary || nnz(active)<Nc
        error('SeMIS:SurvivorSeeds','Distinct survivor selection requires enough binary-weight candidates.');
    end
    candidates=find(active);
    idx=candidates(randperm(numel(candidates),Nc));
    actualSelection='survivors'; candidateCount=numel(candidates);
elseif ~fixed && strcmp(selection,'auto') && all(w==0 | w==1) && nnz(w)>=Nc
    candidates = find(w>0);
    idx = candidates(randperm(numel(candidates),Nc));
    actualSelection='survivors'; candidateCount=numel(candidates);
else
    idx = randsample(numel(w),Nc,true,w);
    actualSelection='weighted'; candidateCount=nnz(w>0);
    candidates=find(w>0);
end
% Preserve the actual draw and its integration-space parents before any
% conditional lift. Recording indices does not draw or reorder samples.
seedInfo=struct('parentIndices',reshape(idx,1,[]),'parentN',size(x,2), ...
    'selection',actualSelection,'requestedSelection',selection, ...
    'candidateCount',candidateCount,'candidateIndices',reshape(candidates,1,[]), ...
    'integrationX',x(:,idx),'integrationY',y(:,idx));
x_sed = x(:,idx);
space=ISDSpace(intBay);
seedInfo.integrationPreserved=strcmp(space.Kind,'identity') && ...
    isempty(space.Lift) && isempty(space.Project);
if ismethod(intBay,'SeedTransportInvariant')
    invariant=intBay.SeedTransportInvariant();
    assert(islogical(invariant) && isscalar(invariant), ...
        'SeMIS:SeedTransportInvariant','SeedTransportInvariant must return a logical scalar.');
    seedInfo.integrationPreserved=seedInfo.integrationPreserved || invariant;
end
if ~isempty(space.Lift) || ~isempty(space.Project) || ~strcmp(space.Kind,'identity')
    selectedG = intBay.G;
    [x_sed,y_sed,intBay,Ncal] = space.LiftSamples(intBay,x_sed,y(:,idx));
    if ~isequaln(selectedG,intBay.G)
        error('SeMIS:LiftSchedule','Lifting seeds must not change the selected ISD parameters.');
    end
    validateattributes(Ncal,{'numeric'},{'scalar','integer','nonnegative','finite'});
    if ~isequal(size(x_sed),[intBay.Ndim,Nc]) || ...
            ~isequal(size(y_sed),[intBay.Nfun,Nc])
        error('SeMIS:LiftLayout','LiftSamples must return one sampling-space seed per selected integration point.');
    end
else
    y_sed = UpdY(intBay,x_sed,y(:,idx));
end
end
