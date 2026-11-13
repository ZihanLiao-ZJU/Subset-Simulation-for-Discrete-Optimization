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
SamplerClass = @aCS;         % CW_aMH, aCS, HMC or ES
simulation = SeMCS(SamplerClass(model));
simulation.NumSam = 512;     % Constant samples per level

%% Run optimization and extract the best designs
raw = simulation.RunIte();
[solution,metric,constraints] = raw.intBay.paretoArchive();
archive_solution = raw.intBay.CacheZ;
archive_objectives = raw.intBay.CacheF;
archive_constraints = raw.intBay.CacheC;

%% Results
result = struct('problem',problem.Name,'solution',solution,'metric',metric, ...
    'constraints',constraints,'sampler',func2str(SamplerClass), ...
    'evaluations',sum(raw.Ncal),'iterations',raw.Nite, ...
    'converged',logical(raw.flg_cvg),'stop_reason',raw.stop_reason, ...
    'archive_solution',archive_solution,'archive_objectives',archive_objectives, ...
    'archive_constraints',archive_constraints,'raw',raw);
