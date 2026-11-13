%% Problem definition
problem = QAPBenchmarks('nug20');
% problem = DiscreteBenchmarks('onemax');
% problem = DiscreteBenchmarks('grid_zdt1','Dimension',20,'Levels',21);

%% Optimization model
model = SuS_Discrete_Optim(problem,problem.Xmin,problem.Xmax, ...
    'StopRule','stability');
% For a multiobjective problem, replace the model above with:
% model = SuS_Multi_Optim(problem,problem.Xmin,problem.Xmax);

%% Sampler and simulation settings
SamplerClass = @aCS;         % aCS, CW_aMH, CW_SMS, ES_2pi, ES_pi or ES_SourceV3
model.Nsam = 512;            % Constant samples per level
model.G.p = 0.1;             % Conditional probability
sampler = SamplerClass();
sampler.BayStc = model;
simulation = SeMCS(sampler);
simulation.MaxIte = 10000;
simulation.Verbose = false;

%% Run optimization
raw = simulation.RunIte();
archive_solution = raw.intBay.CacheZ;
archive_objectives = raw.intBay.CacheF;

%% Extract the best scalar design
[metric,index] = min(archive_objectives);
solution = archive_solution(:,index);
archive_constraints = zeros(0,size(archive_solution,2));
constraints = zeros(0,size(solution,2));

% For a multiobjective problem, replace the four lines above with:
% archive_constraints = raw.intBay.CacheC;
% [solution,metric,constraints] = raw.intBay.paretoArchive();

%% Results
result = struct('problem',problem.Name,'solution',solution,'metric',metric, ...
    'constraints',constraints,'sampler',func2str(SamplerClass), ...
    'evaluations',sum(raw.Ncal),'iterations',raw.Nite, ...
    'converged',logical(raw.flg_cvg),'stop_reason',raw.stop_reason, ...
    'archive_solution',archive_solution,'archive_objectives',archive_objectives, ...
    'archive_constraints',archive_constraints,'raw',raw);
