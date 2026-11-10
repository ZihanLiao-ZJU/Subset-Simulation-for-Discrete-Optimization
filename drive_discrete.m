function result=drive_discrete(example,varargin)
%DRIVE_DISCRETE Discrete optimization with the shared SeMCS and MCMC modules.
% result=drive_discrete;  % QAPLIB nug20: 20 variables, 20 choices each.
% result=drive_discrete('grid_zdt1','Dimension',20,'Levels',21,'MaxIte',30);
% Classical small cases: onemax, knapsack, pmedian, tsp, lotz, cocz, biknapsack.
% Integer-coded grids: grid_rastrigin, grid_rosenbrock, grid_ackley,
% grid_zdt1, grid_dtlz2. QAPLIB: nug20, nug30, sko49, sko100a.
% A constructed objective provider can also be supplied as the first input.
% MaxIte=[] keeps the shared engine default. Patience and MinimumLevels are
% optional method settings; empty values retain the ISD's own defaults.
% Reference oracles are never used during optimization.
if nargin<1,example='nug20';end
activateToolbox(fileparts(mfilename('fullpath')));
p=inputParser;
p.addParameter('NumSam',512,@(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=2)&&all(x==fix(x)));
p.addParameter('MaxIte',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=1&&x==fix(x)));
p.addParameter('ConPrb',.1,@(x)isnumeric(x)&&isscalar(x)&&isreal(x)&&x>0&&x<1);
p.addParameter('Seed',1,@(x)isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=0&&x==fix(x));
p.addParameter('Sampler','aCS',@(x)ischar(x)||isstring(x));
p.addParameter('Verbose',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('Dimension',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=1&&x==fix(x)));
p.addParameter('Levels',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=2&&x==fix(x)));
p.addParameter('Patience',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=1&&x==fix(x)));
p.addParameter('MinimumLevels',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=1&&x==fix(x)));
p.parse(varargin{:});opt=p.Results;
if ischar(example)||(isstring(example)&&isscalar(example))
    if any(strcmpi(char(example),{'nug20','nug30','sko49','sko100a'}))
        if ~isempty(opt.Dimension)||~isempty(opt.Levels)
            error('drive_discrete:FixedInstance','QAP instances prescribe their dimension and integer options.');
        end
        problem=QAPBenchmarks(example);
    else
        problem=DiscreteBenchmarks(example,'Dimension',opt.Dimension,'Levels',opt.Levels);
    end
else
    if ~isempty(opt.Dimension)||~isempty(opt.Levels)
        error('drive_discrete:Provider','Set dimension and integer bounds when constructing the provider.');
    end
    problem=example;
end
options={'Unique',problem.Unique,'PermutationInvariant',problem.PermutationInvariant,'Mapping',problem.Mapping};
for key={'Patience','MinimumLevels'}
    name=key{1};
    if ~isempty(opt.(name)),options=[options,{name,opt.(name)}];end %#ok<AGROW>
end
if problem.Nfun==1
    model=SuS_Discrete_Optim(problem,problem.Xmin,problem.Xmax,options{:},'StopRule','stability');
else
    model=SuS_Multi_opt(problem,problem.Xmin,problem.Xmax,options{:});
end
model.Nsam=opt.NumSam;model.G.p=opt.ConPrb;
kernel=validatestring(char(opt.Sampler),{'aCS','CW_aMH','CW_SMS','ES_pi','ES_2pi','ES_SourceV3'});
sampler=feval(kernel);sampler.BayStc=model;
simulation=SeMCS(sampler);
if ~isempty(opt.MaxIte),simulation.MaxIte=opt.MaxIte;end
simulation.Verbose=opt.Verbose;
rng(opt.Seed,'twister');raw=simulation.RunIte();
Z=raw.intBay.CacheZ;F=raw.intBay.CacheF;C=zeros(0,size(Z,2));
if problem.Ncst,C=raw.intBay.CacheC;end
if problem.Nfun==1
    [metric,index]=min(F);solution=Z(:,index);constraint=C(:,index);
else
    [solution,metric,constraint]=raw.intBay.paretoArchive();
end
result=struct('problem',problem.Name,'solution',solution,'metric',metric, ...
    'constraints',constraint,'sampler',kernel,'seed',opt.Seed, ...
    'evaluations',sum(raw.Ncal),'iterations',raw.Nite, ...
    'converged',logical(raw.flg_cvg),'stop_reason',raw.stop_reason, ...
    'archive_solution',Z,'archive_objectives',F,'archive_constraints',C,'raw',raw);
end

function activateToolbox(root)
% Select this distribution's identical shared engine and support definitions.
folders={'MCMC-samplers','Objective-functions','Support-Fun','Operations on Logarithmic values'};
addpath(root,'-begin');
definitions=dir(fullfile(root,'*.m'));
for k=1:numel(folders)
    folder=fullfile(root,folders{k});addpath(folder,'-begin');
    definitions=[definitions;dir(fullfile(folder,'*.m'))]; %#ok<AGROW>
end
loaded=inmem('-completenames');
for k=1:numel(definitions)
    expected=fullfile(definitions(k).folder,definitions(k).name);
    [~,name]=fileparts(expected);resolved=which(name);
    if ~sameFile(expected,resolved)
        error('DiscreteSuS:PathShadowed','%s resolves outside this toolbox. Change to the toolbox folder and start a fresh MATLAB session.',name);
    end
    for j=1:numel(loaded)
        [~,loadedName,extension]=fileparts(loaded{j});
        if strcmp(name,loadedName)&&strcmpi(extension,'.m')&&~sameFile(expected,loaded{j})
            error('DiscreteSuS:LoadedDefinition','%s was loaded from another toolbox. Start a fresh MATLAB session.',name);
        end
    end
end
end
function yes=sameFile(a,b)
yes=strcmpi(strrep(a,'\','/'),strrep(b,'\','/'));
end
