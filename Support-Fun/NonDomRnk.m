function ranks = NonDomRnk(F, constraints)
%Nondominated ranks for minimization (first front is zero).
% F is objectives-by-points. Optional constraint rows are feasible at <=0.
% Feasible points dominate infeasible points; infeasible points are ordered
% by summed positive violation. Equal vectors have the same rank.
validateattributes(F, {'numeric'}, {'2d','real','finite'});
n = size(F, 2);
if nargin < 2 || isempty(constraints)
    violation = zeros(1, n);
else
    validateattributes(constraints, {'numeric'}, {'2d','real','finite'});
    if isvector(constraints) && numel(constraints) == n
        constraints = reshape(constraints, 1, n);
    end
    if size(constraints, 2) ~= n
        error('DiscreteSuS:ConstraintShape', 'Constraints need one column per objective point.');
    end
    violation = sum(max(constraints, 0), 1);
    if any(~isfinite(violation))
        error('DiscreteSuS:ConstraintOverflow', 'Total constraint violation must be finite.');
    end
end
ranks = zeros(1, n);
if n == 0, return; end
if size(F, 1) < 1
    error('DiscreteSuS:ObjectiveShape', 'At least one objective is required.');
end
dominates = false(n, n);
for i = 1:n
    if violation(i) == 0
        better = all(F(:, i) <= F, 1) & any(F(:, i) < F, 1);
        dominates(i, :) = violation > 0 | (violation == 0 & better);
    else
        dominates(i, :) = violation > violation(i);
    end
end
count = sum(dominates, 1);
assigned = false(1, n);
front = find(count == 0);
level = 0;
while ~isempty(front)
    ranks(front) = level;
    assigned(front) = true;
    count = count - sum(dominates(front, :), 1);
    front = find(~assigned & count == 0);
    level = level + 1;
end
assert(all(assigned), 'DiscreteSuS:Dominance', 'The dominance relation must be acyclic.');
end
