function [logC,yNorm,logH,info] = NormalizeMIS(intBay,x,y,g,logC0)
% Estimate normalizers from RAW per-level [logLi; logpi; model payload].
% logC0 is the externally supplied initial log normalizer (default C_0=1).
% Only proposals actually sampled enter the MIS mixture: the unused next
% proposal g{Nite+1} is irrelevant. The initial factor is shared by all levels.
% An ISD may provide KnownLogNormalizer() for its current G. A finite value
% is the exact ABSOLUTE log integral of h=Li*pi; [] requests estimation.
% intBay may be a cell of frozen integration views, one per sampled level.
% Exact anchors reset accumulated normalizer noise. info.jacobian records
% d log C_i / d log H_m for the common CoV calculation.
if nargin<5, logC0 = 0; end
assert(isnumeric(logC0) && isscalar(logC0) && isreal(logC0) && isfinite(logC0), ...
    'SeMIS:InitialNormalizer','The initial log normalizer must be a finite real scalar.');
Nite = numel(y);
assert(Nite>0 && numel(x)==Nite && numel(g)>=Nite, ...
    'SeMIS:MISSize','Invalid sample/proposal history.');
logC = zeros(Nite,1);
logC(1) = logC0;
known = false(Nite,1); exact = nan(Nite,1);
for k=1:Nite
    probe=MISViewAt(intBay,g,k);
    if ismethod(probe,'KnownLogNormalizer')
        value=probe.KnownLogNormalizer();
        assert(isnumeric(value) && isreal(value) && ...
            (isempty(value) || (isscalar(value) && isfinite(value))), ...
            'SeMIS:KnownNormalizer','KnownLogNormalizer must return [] or a finite real scalar.');
        if isempty(value), continue; end
        known(k)=true; exact(k)=value;
    end
end
if known(1)
    assert(abs(exact(1)-logC0)<=64*eps(max([1,abs(exact(1)),abs(logC0)])), ...
        'SeMIS:KnownNormalizer','Exact first-layer normalizer disagrees with initial metadata.');
    logC(1)=exact(1);
end
B=zeros(Nite,max(Nite-1,0));
yNorm = y;
logH = cell(max(Nite-1,0),1);
for k = 1:Nite-1
    nextView=MISViewAt(intBay,g,k+1);
    [~,~,logH{k}] = ResRat(nextView,x{k},y{k},g{k+1});
    if known(k+1)
        logC(k+1)=exact(k+1);
    else
        increment = logmean(logH{k},'all');
        assert(isreal(increment) && isfinite(increment), ...
            'SeMIS:Normalizer','No finite positive overlap with the next sampled proposal.');
        logC(k+1) = logC(k) + increment;
        B(k+1,:)=B(k,:); B(k+1,k)=1;
    end
end
for k = 1:Nite
    yNorm{k}(1,:) = y{k}(1,:) - logC(k);
end
info=struct('known',known,'exactLogC',exact,'jacobian',[]);
% [] keeps the original CoV arithmetic and default trajectories unchanged.
if any(known), info.jacobian=B; end
end
