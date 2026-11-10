function [u,y,Ncal,baystc] = SamGen(baystc,N)
% Initial ISD: independently sampled normalized Gaussian pi_0 and Li_0=1.
validateattributes(N,{'numeric'},{'scalar','integer','positive','finite'});
[mu,L] = GaussianPrior(baystc);
u = mu + L*randn(baystc.Ndim,N);
baystc = baystc.EvlY(u);
y = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
Ncal = baystc.Ncal;
if size(y,1) ~= baystc.Nfun || size(y,2) ~= N
    error('SeMIS:SampleLayout','Nfun must equal 2 + size(EvlY payload,1).');
end
if any(~isfinite(y(1,:))) || any(abs(y(1,:))>1e-10)
    error('SeMIS:InitialISD','The initial ISD must have log intermediate likelihood zero. Add a reference level before adaptive levels.');
end
z = L\(u-mu);
logpi = -0.5*(baystc.Ndim*log(2*pi)+sum(z.^2,1))-sum(log(diag(L)));
if any(~isfinite(y(2,:))) || any(abs(y(2,:)-logpi)>1e-8*(1+abs(logpi)))
    error('SeMIS:PriorMismatch','EvlPDF must match the Gaussian reference returned by GaussianPrior.');
end
end
