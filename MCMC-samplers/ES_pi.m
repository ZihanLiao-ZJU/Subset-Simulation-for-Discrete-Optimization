classdef ES_pi < ES_2pi
    % Legacy entry point. Uses the full elliptical slice update inherited
    % from ES_2pi; coordinate-wise half ellipses are not valid for a general
    % correlated Gaussian reference. Kept so existing driver choices work.
end
