function [lower,upper] = InitialNormalizerCov(conditionalLower,conditionalUpper,initial)
% Propagate one independent initial normalizer as a COMMON multiplicative
% factor. For independent A,B, CV(AB)^2=CV(A)^2+(1+CV(A)^2)*CV(B)^2.
% Applying the factor once after MisCov prevents artificial dilution across
% levels. The supplied CoVs remain plug-in bounds under MisCov's assumptions.
initial = MISInitialState(struct('initial',initial));
if initial.cov==0
    lower = conditionalLower;
    upper = conditionalUpper;
    return
end
scale = hypot(1,initial.cov);
lower = hypot(initial.cov,scale.*conditionalLower);
upper = hypot(initial.cov,scale.*conditionalUpper);
end
