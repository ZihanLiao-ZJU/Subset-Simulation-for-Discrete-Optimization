classdef aCS
    % Adaptive conditional sampling (pCN-type MCMC) with Nataf transformation
    properties(Constant)
        % SedSav  : use provided seeds at k=1
        % AcpOpt  : target acceptance rate
        % NumIte  : number of adaptive rounds
        % Lambda  : initial proposal scaling
        SedSav = false;
        AcpOpt = 0.44;
        NumIte = 5;
        Lambda = 0.6;
    end

    properties
        % BayStc : Bayesian structure (likelihood, prior, counters, etc.)
        BayStc
        % Default unit base scale and carried lambda; source policies may opt out.
        InitialCorrelation = []  % initial rho; a carried scale takes precedence
        InitialScale = []        % optional initial lambda; carried state takes precedence
        AdaptChains = []         % optional number of chains per adaptive group
        AdaptProposals = []      % update after whole chains accumulate this many proposals
        SeedMode = 'move'         % move | retain | reevaluate first chain state
        ScaleMode = 'unit'        % unit | seeds (standard deviations after whitening)
        CarryScale = true        % reuse StateLambda between sampling calls
        AdaptStatistic = 'realized' % realized acceptance | MH probability
        AdaptParameter = 'lambda' % lambda (innovation scale) | rho (correlation)
        AdaptExponent = 0.5      % diminishing gain at update i: i^(-exponent)
        MaxCorrelation = 0.99   % upper bound used only for rho adaptation
        AdaptFraction = []       % group size ceil(value*Nchains/Nretained)
        AdaptPartial = true      % adapt after an incomplete final group
        StateLambda = []         % last returned scale, read only if CarryScale
        % Optional observable-based diagnostics / between-level score update.
        % Factory freezes a probe from cached seed records, without changing
        % the ISD or Gaussian reference. Off preserves the original path.
        MixingMode = 'off'       % off | diagnostic | score
        MixingProbeFactory = []  % @(bay,seedU,seedY) -> Evaluate and Metadata
        MixingGain = 0.1         % maximum angular gradient step (radians)
        MixingMinProposals = 250
        MixingMinJumps = 20      % also the minimum reward effective sample size
        MixingMaxRatio = 1.25    % limit both sigma and rho changes between calls
    end

    methods
        function obj = ResetState(obj)
            obj.StateLambda = [];
        end

        function [u,y,Ncal,baystc,sampler,info] = SamGen(obj,u_sed,y_sed,steps)
            [u,y,Ncal,baystc,sampler,info] = obj.SamGenState(u_sed,y_sed,steps);
        end

        function [u,y,Ncal,baystc,sampler,info] = SamGenState(obj,u_sed,y_sed,steps)
            % Stateful protocol. Explicit reset-only rectangular sampling
            % retains the historical RNG/arithmetic path below.
            obj.ValidatePolicies;
            Nc = size(u_sed,2);
            validateattributes(steps,{'numeric'},{'real','vector','positive','integer','finite'});
            if Nc<1 || (~isscalar(steps) && numel(steps)~=Nc)
                error('aCS:ChainLengths','Supply a scalar length or one positive length per seed.');
            end
            if size(u_sed,1)~=obj.BayStc.Ndim || ...
                    ~isequal(size(y_sed),[obj.BayStc.Nfun,Nc])
                error('aCS:Seeds','Seed coordinates and payload have incompatible shapes.');
            end
            legacy = isscalar(steps) && strcmp(obj.SeedMode,'move') && ...
                strcmp(obj.ScaleMode,'unit') && ~obj.CarryScale && ...
                strcmp(obj.AdaptStatistic,'realized') && ...
                isempty(obj.AdaptFraction) && obj.AdaptPartial && ...
                isempty(obj.InitialScale) && isempty(obj.AdaptProposals) && ...
                strcmp(obj.AdaptParameter,'lambda') && obj.AdaptExponent==0.5 && ...
                strcmp(obj.MixingMode,'off');
            if legacy
                [u,y,Ncal,baystc,sampler,info] = obj.LegacySamGen(u_sed,y_sed,steps);
            else
                [u,y,Ncal,baystc,sampler,info] = obj.PolicySamGen(u_sed,y_sed,steps);
            end
        end
    end

    methods (Access=private)
        function [u,y,Ncal,baystc,sampler,info] = LegacySamGen(obj,u_sed,y_sed,Ns)
            % Generate MCMC samples targeting prior × likelihood
            % pCN proposal + adaptive scaling on acceptance rate

            % initialization
            baystc = obj.BayStc;
            sedsav = obj.SedSav;
            Ndim = obj.BayStc.Ndim;
            Nfun = obj.BayStc.Nfun;
            Nite = obj.NumIte;
            lambda = obj.Lambda;
            acpopt = obj.AcpOpt;
            firstRho = obj.InitialCorrelation;
            groupSize = obj.AdaptChains;
            if ~isnumeric(firstRho) || ~isreal(firstRho) || ...
                    (~isempty(firstRho) && (~isscalar(firstRho) || ...
                    ~isfinite(firstRho) || firstRho<0 || firstRho>1))
                error('aCS:InitialCorrelation','InitialCorrelation must be empty or a finite scalar in [0,1].');
            end
            if ~isnumeric(groupSize) || ~isreal(groupSize) || ...
                    (~isempty(groupSize) && (~isscalar(groupSize) || ...
                    ~isfinite(groupSize) || groupSize<1 || groupSize~=fix(groupSize)))
                error('aCS:AdaptChains','AdaptChains must be empty or a positive integer.');
            end

            % linear map for correlated Gaussian proposal
            [mu,L] = GaussianPrior(baystc);

            Nc = size(u_sed,2);
            Ncal = 0;

            % allocate memory
            u = zeros(Ndim,Nc,Ns);
            y = zeros(Nfun,Nc,Ns);

            % Gaussian innovations and MH uniforms
            du = pagemtimes(L,randn(Ndim,Nc,Ns));
            p  = rand(1,Nc,Ns);

            % group chains for adaptive updates
            if isempty(groupSize)
                Ngrp = max(ceil(Nc/Nite),1);
            else
                Ngrp = groupSize;
            end
            Nite = ceil(Nc/Ngrp);
            ind_grp = min(repmat([1,Ngrp],Nite,1)+(0:Ngrp:(Nite-1)*Ngrp)',Nc);

            % base proposal scale
            % sigma0 = max(std(u_sed,0,2),eps);
            sigma0 = 1;
            groupAcceptance = zeros(1,Nite);

            for iter = 1:Nite
                % pCN parameters
                sigma = min(lambda*sigma0,1);
                rho   = sqrt(1-sigma.^2);
                if iter==1 && ~isempty(firstRho)
                    % Only the initial proposal changes. Subsequent groups
                    % use the same Lambda-based acceptance-rate adaptation.
                    rho = firstRho;
                    sigma = sqrt(1-rho.^2);
                end
                Nacp  = 0;
                dNcal = 0;
                ind_sam = ind_grp(iter,1):ind_grp(iter,2);

                for k=1:Ns
                    if k == 1 && sedsav
                        % start from seeds
                        u(:,ind_sam,k) = u_sed(:,ind_sam);
                        y(:,ind_sam,k) = y_sed(:,ind_sam);
                    else
                        % pCN proposal
                        u(:,ind_sam,k) = mu + rho.*(u_sed(:,ind_sam)-mu) + sigma.*du(:,ind_sam,k);

                        % evaluate target
                        baystc = baystc.EvlY(u(:,ind_sam,k));
                        y(:,ind_sam,k) = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];

                        % MH accept/reject (log-scale)
                        logratio = y(1,ind_sam,k)-y_sed(1,ind_sam);
                        ind_rej = isnan(logratio) | log(p(1,ind_sam,k)) > min(0,logratio);
                        ind_sam_rej = ind_sam(ind_rej);

                        % revert rejected samples
                        u(:,ind_sam_rej,k) = u_sed(:,ind_sam_rej);
                        y(:,ind_sam_rej,k) = y_sed(:,ind_sam_rej);

                        % statistics
                        Nacp  = Nacp + length(ind_sam)-length(ind_sam_rej);
                        dNcal = dNcal + baystc.Ncal;

                        % update seeds
                        u_sed(:,ind_sam) = u(:,ind_sam,k);
                        y_sed(:,ind_sam) = y(:,ind_sam,k);
                    end
                end

                % acceptance rate
                acp = Nacp/length(ind_sam)/Ns;
                groupAcceptance(iter) = acp;

                % adaptive scaling update
                lambda = exp(log(lambda)+(acp-acpopt)/sqrt(iter));

                % update total calls
                Ncal = Ncal+dNcal;
            end
            baystc.X = reshape(u,Ndim,[]);
            baystc.Y = reshape(y(3:end,:,:),Nfun-2,[]);
            baystc.Ncal = Ncal;
            sampler = obj;
            sampler.BayStc = baystc;
            sampler.StateLambda = lambda;
            info = struct('ChainLengths',[],'Ngenerated',Nc*Ns, ...
                'Nproposals',Nc*Ns,'NseedEvaluations',0,'Nretained',Nc*Ns, ...
                'proposalCounts',repmat(Ns,1,Nc), ...
                'InitialLambda',obj.Lambda,'FinalLambda',lambda, ...
                'ScaleMode',obj.ScaleMode,'SeedMode',obj.SeedMode, ...
                'GroupSize',Ngrp,'GroupAcceptance',groupAcceptance, ...
                'AdaptationCount',Nite,'Sigma0',sigma0);
        end

        function [u,y,Ncal,baystc,sampler,info] = PolicySamGen(obj,u_sed,y_sed,steps)
            baystc = obj.BayStc;
            Ndim = baystc.Ndim; Nfun = baystc.Nfun;
            Nc = size(u_sed,2); rectangular = isscalar(steps);
            if rectangular, lengths = repmat(steps,1,Nc); else, lengths = reshape(steps,1,[]); end
            Nretained = sum(lengths);
            retain = ~strcmp(obj.SeedMode,'move');
            proposalCounts = lengths-double(retain);
            [mu,L] = GaussianPrior(baystc);
            z_sed = L\(u_sed-mu);
            if strcmp(obj.ScaleMode,'seeds')
                sigma0 = max(std(z_sed,0,2),eps);
            else
                sigma0 = ones(Ndim,1);
            end
            lambda = obj.Lambda;
            if ~isempty(obj.InitialScale), lambda = obj.InitialScale; end
            carriedScale = obj.CarryScale && ~isempty(obj.StateLambda);
            if carriedScale, lambda = obj.StateLambda; end
            initialLambda = lambda;
            mixingEnabled = ~strcmp(obj.MixingMode,'off');
            scoreMode = strcmp(obj.MixingMode,'score');
            if mixingEnabled
                probe = obj.MixingProbeFactory(baystc,u_sed,y_sed);
                if ~isstruct(probe) || ~isscalar(probe) || ...
                        ~all(isfield(probe,{'Evaluate','Metadata'})) || ...
                        ~isa(probe.Evaluate,'function_handle') || ...
                        ~isstruct(probe.Metadata) || ~isscalar(probe.Metadata)
                    error('aCS:MixingProbe','The factory must return Evaluate and scalar Metadata.');
                end
                mixing = struct('Mode',obj.MixingMode,'ProbeMetadata',probe.Metadata, ...
                    'UsedLambda',lambda,'NextLambda',lambda,'ProposalCount',0, ...
                    'PositiveJumps',0,'PositiveJumpChains',0,'JumpSum',0, ...
                    'RewardESS',0,'ScoreSum',0,'RelativeGradient',NaN, ...
                    'CurrentMean',0,'CurrentVariance',NaN,'LagOneEstimate',NaN, ...
                    'Updated',false,'Reason','diagnostic','AngleStep',0, ...
                    'ExpectedStep',0,'RealizedJumpSum',0,'Invalid',false, ...
                    'FixedWithinCall',scoreMode,'StationaryAssumption',false);
                hM2 = 0; jumpSquaredSum = 0;
                if scoreMode && (lambda<sin(.02) || lambda>cos(.02))
                    error('aCS:MixingScale','Score mode requires lambda in [sin(.02),cos(.02)].');
                end
            end
            proposalSchedule = ~isempty(obj.AdaptProposals);
            if proposalSchedule
                % A complete chain is the smallest adaptation unit. This
                % schedule overrides chain-group and fractional schedules.
                groupSize = 1;
            elseif ~isempty(obj.AdaptChains)
                groupSize = obj.AdaptChains;
            elseif ~isempty(obj.AdaptFraction)
                groupSize = max(1,ceil(obj.AdaptFraction*Nc/Nretained));
            else
                groupSize = max(ceil(Nc/obj.NumIte),1);
            end
            nGroups = ceil(Nc/groupSize);
            groupAcceptance = nan(1,nGroups);
            groupAdapted = false(1,nGroups);
            adaptationProposals = zeros(1,nGroups);
            chainAcceptance = nan(1,Nc);
            chainProbability = nan(1,Nc);
            starts = [1,1+cumsum(lengths(1:end-1))];
            uf = zeros(Ndim,Nretained); yf = zeros(Nfun,Nretained);
            Ncal = 0; nSeedEvaluations = 0; nAdaptations = 0;
            pendingProposals = 0; pendingAccepted = 0; pendingProbability = 0;
            for ig = 1:nGroups
                ids = (ig-1)*groupSize+1:min(ig*groupSize,Nc);
                sigma = min(lambda*sigma0,1); rho = sqrt(1-sigma.^2);
                if ~carriedScale && (ig==1 || (proposalSchedule && nAdaptations==0)) && ...
                        ~isempty(obj.InitialCorrelation)
                    rho = obj.InitialCorrelation*ones(Ndim,1);
                    sigma = sqrt(1-rho.^2);
                end
                for ic = ids
                    currentU = u_sed(:,ic); currentZ = z_sed(:,ic); currentY = y_sed(:,ic);
                    if mixingEnabled, positiveChain = false; end
                    firstMove = 1;
                    if retain
                        if strcmp(obj.SeedMode,'reevaluate')
                            baystc = baystc.EvlY(currentU);
                            currentY = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                            Ncal = Ncal+baystc.Ncal;
                            nSeedEvaluations = nSeedEvaluations+1;
                        end
                        uf(:,starts(ic)) = currentU; yf(:,starts(ic)) = currentY;
                        firstMove = 2;
                    end
                    nAccepted = 0; sumProbability = 0;
                    for k = firstMove:lengths(ic)
                        if mixingEnabled
                            innovation = randn(Ndim,1);
                            candidateZ = rho.*currentZ+sigma.*innovation;
                        else
                            candidateZ = rho.*currentZ+sigma.*randn(Ndim,1);
                        end
                        candidateU = mu+L*candidateZ;
                        baystc = baystc.EvlY(candidateU);
                        candidateY = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                        Ncal = Ncal+baystc.Ncal;
                        logratio = candidateY(1)-currentY(1);
                        logalpha = min(0,logratio);
                        if isnan(logratio), logalpha = -inf; end
                        alpha = exp(logalpha);
                        accept = log(rand)<=logalpha;
                        sumProbability = sumProbability+alpha;
                        if mixingEnabled
                            % Compute before changing the current state. MH
                            % alpha, rather than a raw crossing, includes
                            % rejection under soft likelihoods as well.
                            currentH = obj.ProbeValue(probe,currentU,currentY);
                            mixing.ProposalCount = mixing.ProposalCount+1;
                            if isfinite(currentH)
                                delta = currentH-mixing.CurrentMean;
                                mixing.CurrentMean = mixing.CurrentMean+delta/mixing.ProposalCount;
                                hM2 = hM2+delta*(currentH-mixing.CurrentMean);
                            else
                                mixing.Invalid = true;
                            end
                            reward = 0; squaredJump = 0;
                            if alpha>0
                                candidateH = obj.ProbeValue(probe,candidateU,candidateY);
                                squaredJump = (candidateH-currentH)^2;
                                reward = alpha*squaredJump;
                                if ~isfinite(reward)
                                    mixing.Invalid = true; reward = 0; squaredJump = 0;
                                end
                            end
                            mixing.JumpSum = mixing.JumpSum+reward;
                            jumpSquaredSum = jumpSquaredSum+reward^2;
                            mixing.RealizedJumpSum = mixing.RealizedJumpSum+double(accept)*squaredJump;
                            if reward>0
                                mixing.PositiveJumps = mixing.PositiveJumps+1;
                                positiveChain = true;
                            end
                            if scoreMode
                                % Derivative of log q_theta(z'|z), keeping
                                % z,z' fixed, where lambda=sin(theta). This
                                % is a local proposal derivative, not a
                                % finite-chain / final MIS variance gradient.
                                angularScore = (rho(1)/sigma(1))*(sum(innovation.^2)-Ndim) ...
                                    -currentZ.'*innovation;
                                scoreReward = reward*angularScore;
                                if ~isfinite(scoreReward)
                                    mixing.Invalid = true;
                                else
                                    mixing.ScoreSum = mixing.ScoreSum+scoreReward;
                                end
                            end
                        end
                        if accept
                            currentU = candidateU; currentZ = candidateZ; currentY = candidateY;
                            nAccepted = nAccepted+1;
                        end
                        idx = starts(ic)+k-1;
                        uf(:,idx) = currentU; yf(:,idx) = currentY;
                    end
                    if mixingEnabled
                        mixing.PositiveJumpChains = mixing.PositiveJumpChains+double(positiveChain);
                    end
                    if proposalCounts(ic)>0
                        chainAcceptance(ic) = nAccepted/proposalCounts(ic);
                        chainProbability(ic) = sumProbability/proposalCounts(ic);
                    end
                    if proposalSchedule
                        pendingProposals = pendingProposals+proposalCounts(ic);
                        pendingAccepted = pendingAccepted+nAccepted;
                        pendingProbability = pendingProbability+sumProbability;
                    end
                end
                if proposalSchedule
                    eligible = pendingProposals>0 && ...
                        (pendingProposals>=obj.AdaptProposals || (obj.AdaptPartial && ig==nGroups));
                    if eligible
                        if strcmp(obj.AdaptStatistic,'probability')
                            groupAcceptance(ig) = pendingProbability/pendingProposals;
                        else
                            groupAcceptance(ig) = pendingAccepted/pendingProposals;
                        end
                        adaptationProposals(ig) = pendingProposals;
                    end
                else
                    if strcmp(obj.AdaptStatistic,'probability')
                        chainStat = chainProbability(ids);
                    else
                        chainStat = chainAcceptance(ids);
                    end
                    valid = isfinite(chainStat);
                    if any(valid), groupAcceptance(ig) = mean(chainStat(valid)); end
                    % No-proposal chains provide no acceptance observation.
                    complete = numel(ids)==groupSize;
                    eligible = any(valid) && (obj.AdaptPartial || complete);
                    if ~obj.AdaptPartial && retain && proposalCounts(ids(end))==0
                        eligible = false;
                    end
                    if eligible, adaptationProposals(ig) = sum(proposalCounts(ids)); end
                end
                if eligible && ~scoreMode
                    nAdaptations = nAdaptations+1;
                    if strcmp(obj.AdaptParameter,'rho')
                        % A scalar correlation is meaningful here because
                        % rho adaptation requires an isotropic unit scale.
                        % Keep lambda as the common state carried to the next
                        % level; the proposal and its MH density are unchanged.
                        gain=nAdaptations^(-obj.AdaptExponent);
                        nextRho=min(rho(1)*exp(-gain*(groupAcceptance(ig)-obj.AcpOpt)), ...
                            obj.MaxCorrelation);
                        lambda=sqrt(1-nextRho^2);
                    elseif obj.AdaptExponent==0.5
                        % Preserve the established default arithmetic path.
                        lambda = exp(log(lambda)+(groupAcceptance(ig)-obj.AcpOpt)/sqrt(nAdaptations));
                    else
                        lambda=exp(log(lambda)+(groupAcceptance(ig)-obj.AcpOpt)* ...
                            nAdaptations^(-obj.AdaptExponent));
                    end
                    groupAdapted(ig) = true;
                    if proposalSchedule
                        % Overshoot belongs to this update; no remainder is
                        % carried across adaptation intervals (Betz Alg.4).
                        pendingProposals = 0; pendingAccepted = 0; pendingProbability = 0;
                    end
                end
            end
            if mixingEnabled
                if mixing.ProposalCount>0
                    mixing.CurrentVariance = hM2/mixing.ProposalCount;
                    mixing.ExpectedStep = mixing.JumpSum/mixing.ProposalCount;
                    if mixing.CurrentVariance>0
                        % Under stationarity this is rho_1(h). With selected
                        % seeds it is a descriptive one-step proxy only.
                        mixing.LagOneEstimate = 1-mixing.ExpectedStep/(2*mixing.CurrentVariance);
                    end
                end
                if jumpSquaredSum>0
                    mixing.RewardESS = mixing.JumpSum^2/jumpSquaredSum;
                end
                if scoreMode && mixing.JumpSum>0 && isfinite(mixing.ScoreSum)
                    mixing.RelativeGradient = mixing.ScoreSum/mixing.JumpSum;
                end
                if scoreMode
                    [lambda,mixing] = obj.UpdateMixingScale(lambda,mixing);
                else
                    mixing.ScoreSum = NaN;
                    mixing.NextLambda = lambda;
                end
            end
            if rectangular
                % Flat buffers contain whole chains consecutively; preserve
                % the established rectangular [dimension,chain,step] API.
                u = permute(reshape(uf,Ndim,steps,Nc),[1,3,2]);
                y = permute(reshape(yf,Nfun,steps,Nc),[1,3,2]);
                chainLengths = [];
            else
                u = uf; y = yf; chainLengths = lengths;
            end
            baystc.X = reshape(u,Ndim,[]);
            baystc.Y = reshape(y(3:end,:,:),Nfun-2,[]);
            baystc.Ncal = Ncal;
            sampler = obj; sampler.BayStc = baystc; sampler.StateLambda = lambda;
            info = struct('ChainLengths',chainLengths, ...
                'Ngenerated',sum(proposalCounts),'Nproposals',sum(proposalCounts), ...
                'NseedEvaluations',nSeedEvaluations,'Nretained',Nretained, ...
                'proposalCounts',proposalCounts,'InitialLambda',initialLambda, ...
                'FinalLambda',lambda,'ScaleMode',obj.ScaleMode,'SeedMode',obj.SeedMode, ...
                'GroupSize',groupSize,'GroupAcceptance',groupAcceptance, ...
                'GroupAdapted',groupAdapted,'AdaptationCount',nAdaptations, ...
                'AdaptationProposalCounts',adaptationProposals(groupAdapted), ...
                'AdaptationChainEnds',min(find(groupAdapted)*groupSize,Nc), ...
                'AdaptProposals',obj.AdaptProposals, ...
                'ChainAcceptance',chainAcceptance,'ChainProbability',chainProbability, ...
                'Sigma0',sigma0);
            if ~strcmp(obj.AdaptParameter,'lambda') || obj.AdaptExponent~=0.5
                info.AdaptParameter=obj.AdaptParameter;
                info.AdaptExponent=obj.AdaptExponent;
                info.MaxCorrelation=obj.MaxCorrelation;
            end
            if mixingEnabled, info.Mixing = mixing; end
        end

        function value = ProbeValue(~,probe,u,y)
            value = probe.Evaluate(u,y);
            if ~isnumeric(value) && ~islogical(value)
                error('aCS:MixingProbe','The observable must return a real scalar.');
            end
            if ~isreal(value) || ~isscalar(value)
                error('aCS:MixingProbe','The observable must return a real scalar.');
            end
            value = double(value);
            % Nonfinite values invalidate the whole call's update; do not
            % silently delete selected proposals from the objective.
        end

        function [lambda,mixing] = UpdateMixingScale(obj,lambda,mixing)
            if mixing.Invalid || ~isfinite(mixing.JumpSum) || ...
                    ~isfinite(mixing.ScoreSum) || ~isfinite(mixing.CurrentVariance) || ...
                    ~isfinite(mixing.RewardESS)
                mixing.Reason = 'invalid';
            elseif mixing.ProposalCount<obj.MixingMinProposals
                mixing.Reason = 'few_proposals';
            elseif mixing.PositiveJumps<obj.MixingMinJumps
                mixing.Reason = 'few_jumps';
            elseif mixing.RewardESS<obj.MixingMinJumps
                mixing.Reason = 'low_reward_ess';
            elseif mixing.CurrentVariance<=eps(max(1,abs(mixing.CurrentMean)))^2
                mixing.Reason = 'constant_probe';
            elseif obj.MixingGain==0
                mixing.Reason = 'fixed_control';
            elseif obj.MixingMaxRatio==1
                mixing.Reason = 'fixed_ratio';
            else
                mixing.RelativeGradient = mixing.ScoreSum/mixing.JumpSum;
                theta = asin(lambda);
                requestedStep = obj.MixingGain*tanh(mixing.RelativeGradient);
                nextTheta = min(pi/2-.02,max(.02,theta+requestedStep));
                % An angular bound alone permits large relative steps near
                % zero. Bound both innovation and correlation multipliers.
                rho = sqrt(1-lambda^2); ratio = obj.MixingMaxRatio;
                lo = max([.02,asin(lambda/ratio),acos(min(1,rho*ratio))]);
                hi = min([pi/2-.02,asin(min(1,lambda*ratio)),acos(rho/ratio)]);
                nextTheta = min(hi,max(lo,nextTheta));
                mixing.AngleStep = nextTheta-theta;
                if mixing.AngleStep~=0
                    lambda = sin(nextTheta);
                    mixing.Updated = true;
                end
                mixing.Reason = 'score';
            end
            mixing.NextLambda = lambda;
        end

        function ValidatePolicies(obj)
            if ~((ischar(obj.AdaptParameter) && isrow(obj.AdaptParameter)) || ...
                    (isstring(obj.AdaptParameter) && isscalar(obj.AdaptParameter))) || ...
                    ~any(strcmp(obj.AdaptParameter,{'lambda','rho'}))
                error('aCS:AdaptParameter','AdaptParameter must be lambda or rho.');
            end
            validateattributes(obj.AdaptExponent,{'numeric'},{'scalar','real','finite','positive','<=',1});
            validateattributes(obj.MaxCorrelation,{'numeric'},{'scalar','real','finite','>',0,'<',1});
            if strcmp(obj.AdaptParameter,'rho') && ~strcmp(obj.ScaleMode,'unit')
                error('aCS:RhoScaleMode','Correlation adaptation requires unit reference scale.');
            end
            if ~((ischar(obj.MixingMode) && isrow(obj.MixingMode)) || ...
                    (isstring(obj.MixingMode) && isscalar(obj.MixingMode))) || ...
                    ~any(strcmp(obj.MixingMode,{'off','diagnostic','score'}))
                error('aCS:MixingPolicy','MixingMode must be off, diagnostic, or score.');
            end
            if ~isempty(obj.MixingProbeFactory) && ~isa(obj.MixingProbeFactory,'function_handle')
                error('aCS:MixingPolicy','MixingProbeFactory must be empty or a function handle.');
            end
            if ~strcmp(obj.MixingMode,'off') && isempty(obj.MixingProbeFactory)
                error('aCS:MixingPolicy','An enabled mixing policy needs a probe factory.');
            end
            validateattributes(obj.MixingGain,{'numeric'},{'scalar','real','finite','nonnegative'});
            validateattributes(obj.MixingMinProposals,{'numeric'},{'scalar','real','finite','positive','integer'});
            validateattributes(obj.MixingMinJumps,{'numeric'},{'scalar','real','finite','positive','integer'});
            validateattributes(obj.MixingMaxRatio,{'numeric'},{'scalar','real','finite','>=',1});
            if strcmp(obj.MixingMode,'score') && ...
                    (~strcmp(obj.ScaleMode,'unit') || ~obj.CarryScale || ~isempty(obj.InitialCorrelation) || ...
                    ~strcmp(obj.AdaptParameter,'lambda') || obj.AdaptExponent~=0.5)
                error('aCS:MixingPolicy','Score mode requires unit scale, CarryScale, default lambda adaptation, and no InitialCorrelation.');
            end
            if ~((ischar(obj.SeedMode) && isrow(obj.SeedMode)) || ...
                    (isstring(obj.SeedMode) && isscalar(obj.SeedMode))) || ...
                    ~any(strcmp(obj.SeedMode,{'move','retain','reevaluate'}))
                error('aCS:SeedMode','SeedMode must be move, retain, or reevaluate.');
            end
            if ~((ischar(obj.ScaleMode) && isrow(obj.ScaleMode)) || ...
                    (isstring(obj.ScaleMode) && isscalar(obj.ScaleMode))) || ...
                    ~any(strcmp(obj.ScaleMode,{'unit','seeds'}))
                error('aCS:ScaleMode','ScaleMode must be unit or seeds.');
            end
            if ~((ischar(obj.AdaptStatistic) && isrow(obj.AdaptStatistic)) || ...
                    (isstring(obj.AdaptStatistic) && isscalar(obj.AdaptStatistic))) || ...
                    ~any(strcmp(obj.AdaptStatistic,{'realized','probability'}))
                error('aCS:AdaptStatistic','AdaptStatistic must be realized or probability.');
            end
            flags = {obj.CarryScale,obj.AdaptPartial};
            for k = 1:numel(flags)
                f = flags{k};
                if ~((islogical(f) || isnumeric(f)) && isreal(f) && ...
                        isscalar(f) && isfinite(f) && (f==0 || f==1))
                    error('aCS:PolicyFlag','CarryScale and AdaptPartial must be scalar logical flags.');
                end
            end
            if ~isnumeric(obj.StateLambda) || ~isreal(obj.StateLambda)
                error('aCS:StateLambda','StateLambda must be empty or a positive finite scalar.');
            elseif ~isempty(obj.StateLambda)
                validateattributes(obj.StateLambda,{'numeric'},{'scalar','real','finite','positive'});
            end
            firstRho = obj.InitialCorrelation;
            if ~isnumeric(firstRho) || ~isreal(firstRho) || ...
                    (~isempty(firstRho) && (~isscalar(firstRho) || ...
                    ~isfinite(firstRho) || firstRho<0 || firstRho>1))
                error('aCS:InitialCorrelation','InitialCorrelation must be empty or a finite scalar in [0,1].');
            end
            groupSize = obj.AdaptChains;
            if ~isnumeric(groupSize) || ~isreal(groupSize) || ...
                    (~isempty(groupSize) && (~isscalar(groupSize) || ...
                    ~isfinite(groupSize) || groupSize<1 || groupSize~=fix(groupSize)))
                error('aCS:AdaptChains','AdaptChains must be empty or a positive integer.');
            end
            if ~isnumeric(obj.InitialScale) || ~isreal(obj.InitialScale) || ...
                    (~isempty(obj.InitialScale) && (~isscalar(obj.InitialScale) || ...
                    ~isfinite(obj.InitialScale) || obj.InitialScale<=0))
                error('aCS:InitialScale','InitialScale must be empty or a positive finite scalar.');
            end
            count = obj.AdaptProposals;
            if ~isnumeric(count) || ~isreal(count) || ...
                    (~isempty(count) && (~isscalar(count) || ~isfinite(count) || ...
                    count<1 || count~=fix(count)))
                error('aCS:AdaptProposals','AdaptProposals must be empty or a positive integer.');
            end
            if ~isnumeric(obj.AdaptFraction) || ~isreal(obj.AdaptFraction)
                error('aCS:AdaptFraction','AdaptFraction must be empty or a positive finite scalar.');
            elseif ~isempty(obj.AdaptFraction)
                validateattributes(obj.AdaptFraction,{'numeric'},{'scalar','real','finite','positive'});
            end
        end
    end
end
