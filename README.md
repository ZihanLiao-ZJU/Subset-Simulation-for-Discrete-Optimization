# Subset Simulation for Discrete Optimization

MATLAB implementations of single-objective and multiobjective subset simulation for bounded integers, subsets and permutations. Both methods use the same `SeMCS` engine and MCMC samplers as reliability estimation and Bayesian updating in the unified framework.

This toolbox builds on the multiobjective subset simulation approach described in:

> Xiao, D., Liao, Z., & Li, B. (2022). A strategy for optimal sensor placement in modal signal-to-noise ratio-based structural health monitoring. In *2022 8th International Conference on Hydraulic and Civil Engineering: Deep Space Intelligent Development and Utilization Forum (ICHCE)* (pp. 1–8). IEEE. [https://doi.org/10.1109/ICHCE57331.2022.10042699](https://doi.org/10.1109/ICHCE57331.2022.10042699)

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Built-in examples](#built-in-examples)
- [Custom problems](#custom-problems)
- [Outputs](#outputs)
- [Citation](#citation)
- [License](#license)

## Overview

The toolbox solves discrete minimization problems:

$$
\min_{\boldsymbol z\in\mathcal D}\ \boldsymbol f(\boldsymbol z),\qquad \boldsymbol c(\boldsymbol z)\le\boldsymbol0.
$$

Here `z` is a discrete design, `f` contains one or more objectives, and `c` contains optional constraints. The algorithms sample standard-normal coordinates and map them to discrete designs. `SuS_Discrete_Optim` adapts objective thresholds; `SuS_Multi_Optim` uses nondominated ranks against a reference set fixed within each level. Evaluated designs are cached, and the result is the best evaluated scalar design or a feasible nondominated archive.

## Requirements

- MATLAB with Statistics and Machine Learning Toolbox; tested with MATLAB R2026a.
- The included examples and QAPLIB data require no additional downloads.

## Repository structure

```text
.
├── drive_discrete.m                     % choose a problem, configure a run and return the result
├── SeMCS.m                              % the shared sequential sampling engine
├── SuS_Discrete_Optim.m                 % the single-objective ISD
├── SuS_Multi_Optim.m                    % the multiobjective ISD
├── MCMC-samplers/                       % the shared sampling algorithms
├── Objective-functions/                 % classical examples and QAPLIB instances
├── Support-Fun/                         % discrete mappings and common framework support
├── Operations on Logarithmic values/    % logarithmic arithmetic helpers
├── LICENSE
└── README.md
```

## Quick start

Start a fresh MATLAB session in the repository root:

```matlab
single = drive_discrete('nug20');
multi = drive_discrete('grid_zdt1','Dimension',20,'Levels',21);

single.solution
single.metric
multi.solution
multi.metric
```

Choose a sampler and an explicit level budget:

```matlab
r = drive_discrete('grid_dtlz2','Dimension',20,'Levels',21, ...
    'NumSam',512,'MaxIte',30,'ConPrb',0.1,'Seed',1, ...
    'Sampler','CW_aMH');
```

Defaults are `nug20`, 512 samples per level, conditional probability 0.1, seed 1 and `aCS`. Omitting `MaxIte` retains the common engine's 10000-level cap. The native stopping rules usually stop earlier: scalar runs monitor the best evaluated objective; multiobjective runs monitor the distinct feasible nondominated objective vectors. `Patience` defaults to 5 and `MinimumLevels` to 3. These are search stopping rules, not certificates of global optimality or complete Pareto coverage.

Available samplers are `aCS`, `CW_aMH`, `CW_SMS`, `ES_2pi`, `ES_pi` and `ES_SourceV3`. `ES_pi` is a compatibility entry for `ES_2pi`; `ES_SourceV3` preserves the separate coordinate-angle source kernel.

## Built-in examples

| Selector | Problem |
|---|---|
| `nug20`, `nug30`, `sko49`, `sko100a` | QAPLIB quadratic assignment instances, with 20, 30, 49 or 100 permutation coordinates |
| `grid_rastrigin`, `grid_rosenbrock`, `grid_ackley` | Integer-coded grids of classical single-objective functions |
| `grid_zdt1`, `grid_dtlz2` | Integer-coded grids of two-objective ZDT1 and three-objective DTLZ2 |
| `onemax`, `knapsack`, `pmedian`, `tsp` | Small binary, subset and permutation examples |
| `lotz`, `cocz`, `biknapsack` | Small multiobjective examples |

For grid examples, `Dimension` selects the number of variables and `Levels` the number of integer choices per variable. The provider maps these integer codes to the physical grid before evaluating the classical function. These are explicitly restricted-grid variants. Use odd `Levels`; Rosenbrock additionally requires `Levels-1` divisible by four. QAPLIB instances have fixed dimensions. Reference calculations in `Oracle()` are for external assessment and are not used by the optimizer.

QAP data follow [QAPLIB](https://coral.ise.lehigh.edu/data-sets/qaplib/). The grid multiobjective examples use the published [ZDT](https://sop.tik.ee.ethz.ch/publicationListFiles/zdt2000a.pdf) and [DTLZ](https://sop.tik.ee.ethz.ch/publicationListFiles/dtlz2002a.pdf) function definitions; LOTZ and COCZ follow [Laumanns, Thiele and Zitzler (2004)](https://sop.tik.ee.ethz.ch/publicationListFiles/ltz2004a.pdf). The small geometric and knapsack instances are self-contained constructed examples.

## Custom problems

Supply one design per column and return one objective value per column. For example:

```matlab
addpath(genpath(pwd));
objective = @(Z) sum((Z-3).^2,1);
model = SuS_Discrete_Optim(objective,zeros(20,1),6*ones(20,1), ...
    'Unique',false,'PermutationInvariant',false, ...
    'Mapping','uniform','StopRule','stability');
model.Nsam = 512;
sampler = aCS;
sampler.BayStc = model;
simulation = SeMCS(sampler);
simulation.MaxIte = 30;
raw = simulation.RunIte();
[bestValue,index] = min(raw.intBay.CacheF);
bestDesign = raw.intBay.CacheZ(:,index);
```

An object or struct can instead provide `EvlF(Z)`; the original `EvlOBF(Z)` interface is also supported. For `SuS_Multi_Optim`, return an objectives-by-designs matrix and declare `Nfun` as the objective count. A bare function handle can declare `ObjectiveCount` in the constructor. Constrained multiobjective providers additionally define `Ncst` and `EvlCST(Z)`, returning constraints-by-designs values feasible at `<=0`. Scalar constraints are encoded in an explicit finite penalty.

`Unique=true` enforces distinct coordinates. Use `PermutationInvariant=true` only for unordered designs whose objectives are unchanged by sorting and whose bounds are common. Permutations use `PermutationInvariant=false`, `Mapping='rank'` and common bounds `1:dimension`. `Mapping='uniform'` gives equal integer bins; `quantile` retains the original rounded Gaussian-CDF mapping. Direct scalar construction retains `StopRule='legacy'` unless stability is requested.

The cache retains the first evaluation of each design; objectives should be deterministic or use explicitly fixed common random numbers.

## Outputs

| Field | Meaning |
|---|---|
| `solution`, `metric`, `constraints` | Best evaluated scalar design/value, or one design per distinct feasible nondominated objective vector |
| `evaluations` | Newly evaluated discrete designs; cached repeats add no calls |
| `iterations`, `converged`, `stop_reason` | Sampled levels and native termination information |
| `archive_solution`, `archive_objectives`, `archive_constraints` | All uniquely evaluated designs and their cached values |
| `raw` | The shared SeMCS history and final ISD |

For multiobjective problems, objectives and constraints at one new design count as one joint evaluation. The final archive includes evaluated proposals rejected by MCMC. `raw.x` contains Gaussian coordinates; `raw.intBay.u2ints(U)` maps them to discrete designs. The optimization driver does not perform evidence or failure-probability postprocessing.

## Citation

If you use this toolbox in your research, please cite:

> Xiao, D., Liao, Z., & Li, B. (2022). A strategy for optimal sensor placement in modal signal-to-noise ratio-based structural health monitoring. In *2022 8th International Conference on Hydraulic and Civil Engineering: Deep Space Intelligent Development and Utilization Forum (ICHCE)* (pp. 1–8). IEEE. [https://doi.org/10.1109/ICHCE57331.2022.10042699](https://doi.org/10.1109/ICHCE57331.2022.10042699)

```bibtex
@inproceedings{xiao2022strategy,
  author    = {Xiao, Dahai and Liao, Zihan and Li, Binbin},
  title     = {A Strategy for Optimal Sensor Placement in Modal Signal-to-noise Ratio-based Structural Health Monitoring},
  booktitle = {2022 8th International Conference on Hydraulic and Civil Engineering: Deep Space Intelligent Development and Utilization Forum (ICHCE)},
  year      = {2022},
  pages     = {1--8},
  publisher = {IEEE},
  doi       = {10.1109/ICHCE57331.2022.10042699}
}
```

## License

Distributed under the [MIT License](LICENSE). Retained third-party notices are supplied with the corresponding support routines.
