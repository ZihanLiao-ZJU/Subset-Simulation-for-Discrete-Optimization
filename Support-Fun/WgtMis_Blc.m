function [w_mis, log_alpha_full] = WgtMis_Blc(intBay,x,y,g,c)
% Balance weights alpha_j(u)=N_j q_j(u)/sum_l N_l q_l(u).
% q_j is normalized using c(j); y may have raw or normalized first rows.
% The second output contains cross-level balance weights, not log densities.
% A cell intBay contains frozen per-proposal integration evaluators.
Nite = numel(y);
assert(numel(x)==Nite && numel(g)>=Nite && numel(c)>=Nite, ...
    'SeMIS:MISSize','One sample set, proposal and normalizer are required per level.');
Nsam = cellfun(@(v) prod(size(v,2:3)), y);
assert(all(Nsam>0),'SeMIS:EmptyLevel','MIS cannot include an empty sample level.');
assert(all(isfinite(c(1:Nite))), 'SeMIS:Normalizer', ...
    'Sampled proposals must have finite log normalizing constants.');
logN = log(Nsam(:));
w_mis = cell(Nite,1);
log_alpha_full = cell(Nite,1);
% An optional value signature certifies equal unnormalized densities within
% one evaluator class. Repeated fixed proposals are evaluated once; every
% sampled level and its original balance weight remain in the estimator.
views=cell(Nite,1);keys=cell(Nite,1);group=1:Nite;
for j=1:Nite
    views{j}=MISViewAt(intBay,g,j);
    if ismethod(views{j},'DensitySignature')
        keys{j}=views{j}.DensitySignature;
        for k=1:j-1
            if strcmp(class(views{j}),class(views{k}))&&c(j)==c(k)&& ...
                    ~isempty(keys{k})&&isequaln(keys{j},keys{k})
                group(j)=group(k);break;
            end
        end
    end
end
for ite = 1:Nite
    logq = zeros(Nite,Nsam(ite));
    for j = 1:Nite
        if group(j)~=j,logq(j,:)=logq(group(j),:);continue;end
        view = views{j};
        view.X = x{ite}(1:end,:);
        view.Y = y{ite}(3:end,:);
        logq(j,:) = view.EvlLKF + view.EvlPDF - c(j);
    end
    assert(isreal(logq) && ~any(isnan(logq(:)) | logq(:)==inf), ...
        'SeMIS:InvalidDensity','Proposal log densities must be real, finite or -Inf.');
    assert(all(isfinite(logq(ite,:))), 'SeMIS:ProposalSupport', ...
        'A sample has zero density under the proposal recorded for its level.');
    weighted = logq + logN;
    alpha = weighted - logsum(weighted,1);
    w_mis{ite} = reshape(alpha(ite,:),size(y{ite},2:3));
    log_alpha_full{ite} = alpha;
end
end
