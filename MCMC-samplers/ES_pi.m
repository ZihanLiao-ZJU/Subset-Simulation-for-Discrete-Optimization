classdef ES_pi < ES_2pi
    % Historical name for the common-angle slice sampler.
    methods
        function obj = ES_pi(model)
            if nargin < 1, model = []; end
            obj@ES_2pi(model);
        end
    end
end
