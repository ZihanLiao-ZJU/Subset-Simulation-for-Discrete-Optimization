function [W,R,psi,A,gamma,logHmean,status]=MISJointInfluence(context,logEstimate,logW,logAlpha,logH,B)
% Shared first-order MIS observables. No dependence closure is imposed here.
L=numel(logW);N=cellfun(@numel,logW);status='';
if isempty(B),B=[zeros(1,L-1);tril(ones(L-1))];end
assert(isequal(size(B),[L,L-1])&&all(isfinite(B),'all'), ...
    'SeMIS:NormalizerJacobian','Invalid first-order normalizer sensitivity.');
W=cell(L,1);R=cell(L,1);psi=cell(L,1);A=zeros(L,L);
logHmean=zeros(max(0,L-1),1);
for i=1:L
    W{i}=exp(logW{i}(:)-logEstimate);
    assert(size(logAlpha{i},2)==N(i),'SeMIS:MISSize','Balance weights do not match particles.');
    A(i,:)=mean(W{i}.*exp(logAlpha{i}).',1);
    if i<L,logHmean(i)=logmean(logH{i},'all');end
    R{i}=zeros(N(i),1);
end
gamma=B.'*sum(A,1).';gamma=[gamma;0];
for i=1:L
    hasSeed=isfield(context,'SeedInfo')&&numel(context.SeedInfo)>=i+1&&~isempty(context.SeedInfo{i+1});
    needsRatio=i<L&&(gamma(i)~=0||hasSeed);
    if needsRatio
        if ~isfinite(logHmean(i)),status='unavailable_required_overlap';return;end
        R{i}=exp(logH{i}(:)-logHmean(i));
    end
    psi{i}=W{i}+gamma(i)*R{i};
end
end
