function weights = MISTargetWeights(context,y,j)
% log(target/q_i) in the common integration measure, using cached payload.
% EvlTarget(j) defines a target density, never a sampling-space density.
if nargin<3, j=1; end
validateattributes(j,{'numeric'},{'scalar','integer','positive','<=',context.TargetCount});
assert(iscell(y) && numel(y)==context.Nite,'SeMIS:MISSize', ...
    'Target weights require one normalized record per integration level.');
weights=cell(context.Nite,1);
for k=1:context.Nite
    yi=y{k}; view=MISViewAt(context.Views,context.G,k);
    if ismethod(view,'EvlTarget')
        view.X=context.X{k}(1:end,:); view.Y=yi(3:end,:);
        target=view.EvlTarget(j);
        assert(isnumeric(target) && isequal(size(target),[1,prod(size(yi,2:3))]) && ...
            isreal(target) && ~any(isnan(target) | target==inf), ...
            'SeMIS:TargetDensity','EvlTarget must return one finite or -Inf log target per integration point.');
        wi=target-yi(1,:)-yi(2,:);
    elseif strcmp(context.Domain,'bayesian')
        assert(size(yi,1)>=4 && context.Descriptors{k}.likelihoodDiagnostics, ...
            'SeMIS:TargetDensity','A custom integration payload requires EvlTarget(j).');
        % Preserve the legacy Bayesian arithmetic exactly.
        wi=yi(3,:)+yi(4,:)-yi(1,:)-yi(2,:);
        wi(isnan(wi))=-inf;
    else
        assert(size(yi,1)>=3+j && strcmp(context.Descriptors{k}.kind,'identity'), ...
            'SeMIS:TargetDensity','A projected reliability target requires EvlTarget(j).');
        wi=yi(3,:)-yi(1,:)-yi(2,:);
        wi(yi(3+j,:)>=0)=-inf;
    end
    assert(isreal(wi) && ~any(isnan(wi) | wi==inf), ...
        'SeMIS:TargetDensity','Target/proposal ratios must be finite or -Inf on retained samples.');
    weights{k}=reshape(wi,size(yi,2:3));
end
end
