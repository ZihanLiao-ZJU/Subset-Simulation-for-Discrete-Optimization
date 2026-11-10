function [mu,L] = GaussianPrior(intBay)
% Gaussian reference pi_i = N(mu,L*L'). Parameters belong to G, so MIS can
% reconstruct historical ISDs. L_pi/mu_pi are aliases from reliability v1.
if ismethod(intBay,'PriorGaussian')
    [mu,L] = intBay.PriorGaussian();
else
    mu = zeros(intBay.Ndim,1);
    L = eye(intBay.Ndim);
    if isstruct(intBay.G)
        if isfield(intBay.G,'mu'), mu = intBay.G.mu;
        elseif isfield(intBay.G,'mu_pi'), mu = intBay.G.mu_pi; end
        if isfield(intBay.G,'L'), L = intBay.G.L;
        elseif isfield(intBay.G,'L_pi'), L = intBay.G.L_pi; end
    end
end
mu = mu(:);
if ~isequal(size(mu),[intBay.Ndim,1]) || ...
        ~isequal(size(L),[intBay.Ndim,intBay.Ndim]) || ...
        any(~isfinite(mu)) || any(~isfinite(L),'all') || ...
        ~istril(L) || any(diag(L)<=0)
    error('SeMIS:GaussianPrior','Provide a finite mean and a lower Cholesky factor with positive diagonal.');
end
end
