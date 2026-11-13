% Adapted from Zihan Liao's public MATLAB MCMC samplers; repository LICENSE applies.
classdef aCS < MCMCBase
    % sigma0=1; updated lambda is returned in b.G.lambda across levels.
    properties
        Lambda = 0.6
        AcpOpt = 0.44
        NumIte = 5
        Adapt = true
    end

    methods

        function obj = aCS(b)
            if nargin < 1
                b = [];
            end

            obj@MCMCBase(b);
        end

        function [u, y, Ncal, b] = SamGen(obj, varargin)
            [b, u0, y0, Ncal, save, L, Nc, Ns] = obj.prepare(varargin{:});
            d = b.Ndim;
            u = zeros(d, Nc, Ns);
            y = zeros(b.Nfun, Nc, Ns);
            lambda = obj.state(b, 'lambda', obj.Lambda);
            if ~isscalar(lambda) || ~isfinite(lambda) || lambda <= 0
                error('aCS:State', 'G.lambda must be finite and positive.');
            end

            width = max(ceil(Nc / obj.NumIte), 1);
            group = 0;

            for start = 1:width:Nc
                group = group + 1;
                ids = start:min(start + width - 1, Nc);
                sigma = min(lambda, 1);
                rho = sqrt(1 - sigma ^ 2);
                first = 1 + save;
                rates = zeros(1,numel(ids));
                for chain = ids
                    current = u0(:,chain);
                    values = y0(:,chain);
                    accepted = 0;
                    if save
                        u(:,chain,1) = current;
                        y(:,chain,1) = values;
                    end
                    for step = first:Ns
                        proposal = rho .* current + sigma .* (L * randn(d,1));
                        [b, proposedValues, calls] = obj.evaluate(b, proposal);
                        Ncal = Ncal + calls;
                        if ~obj.reject(proposedValues(1), values(1), rand)
                            current = proposal;
                            values = proposedValues;
                            accepted = accepted + 1;
                        end
                        u(:,chain,step) = current;
                        y(:,chain,step) = values;
                    end
                    if Ns >= first
                        rates(chain-start+1) = accepted / (Ns-first+1);
                    end
                end
                if obj.Adapt && Ns >= first
                    lambda = exp(log(lambda) + (mean(rates) - obj.AcpOpt) / sqrt(group));
                end
            end

            b.G.lambda = lambda;
        end

    end
end
