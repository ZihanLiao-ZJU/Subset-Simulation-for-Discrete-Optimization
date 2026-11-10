function [lograt,shift,rawLogRatio] = ResRat(intBay,x,y,g)
% lograt is max-shifted for seed selection (legacy Bayesian v3 convention).
% rawLogRatio, NOT lograt alone, must be used for normalizing constants.
% Only rows 1:2 (logLi/logpi) are interpreted by this shared function.
shape = [size(y,2),size(y,3)];
y = reshape(y,size(y,1),[]);
intBay.G = g;
intBay.X = reshape(x,size(x,1),[]);
intBay.Y = y(3:end,:);
raw = intBay.EvlLKF + intBay.EvlPDF - y(1,:) - y(2,:);
if any(isnan(raw)) || any(raw==inf)
    error('SeMIS:InvalidDensityRatio','ISD ratios must be finite or -Inf on the sampled support.');
end
shift = max(raw,[], 'all');
if shift == -inf
    scaled = -inf(size(raw));
else
    scaled = raw-shift;
end
lograt = reshape(scaled,shape);
rawLogRatio = reshape(raw,shape);
end
