function [var_low,var_up] = MoveBlock(log_w,chainLengths)
% Reliability v1 Appendix B: rows are chains, columns are successive states.
% Chain means retain within-chain correlation; moving blocks of those means
% allow dependence between resampled chains. A single chain uses time blocks.
if ~iscell(log_w), log_w = {log_w}; end
if nargin<2 || isempty(chainLengths), chainLengths=cell(size(log_w)); end
if ~iscell(chainLengths), chainLengths={chainLengths}; end
assert(numel(chainLengths)==numel(log_w),'SeMIS:ChainLayout','Chain layout count differs from levels.');
variances = zeros(numel(log_w),1);
for k = 1:numel(log_w)
    a = log_w{k};
    assert(ismatrix(a) && ~isempty(a) && isreal(a) && ...
        ~any(isnan(a(:)) | a(:)==inf), 'SeMIS:InvalidWeights', ...
        'Weights must be a nonempty [chains, steps] matrix of real log values.');
    shift = max(a,[],'all');
    if shift == -inf, continue; end
    w = exp(a-shift);
    lengths=chainLengths{k};
    if ~isempty(lengths)
        validateattributes(lengths,{'numeric'},{'vector','integer','positive','finite'});
        assert(sum(lengths)==numel(w),'SeMIS:ChainLayout','Chain lengths do not sum to retained records.');
        nc=numel(lengths);
        if nc==1
            v=blockVariance(w(:));
        elseif all(lengths==lengths(1))
            % Ragged-format storage is chain-major, unlike the rectangular
            % [chains,steps] layout. Equal lengths recover the original rule.
            v=blockVariance(mean(reshape(w(:),lengths(1),nc),1).');
        else
            % Linearize the ratio sum(chain totals)/sum(chain lengths).
            % Center each chain by its own length so unequal lengths alone
            % do not create spurious variance for a constant integrand.
            ends=[0;cumsum(lengths(:))]; totals=zeros(nc,1);
            for j=1:nc, totals(j)=sum(w(ends(j)+1:ends(j+1))); end
            C=nc/numel(w)*(totals-lengths(:)*mean(w(:)));
            v=blockVariance(C);
        end
    elseif size(w,1)>1
        C = mean(w,2);
        v = blockVariance(C);
    else
        % With one chain, between-chain variance is unavailable: use time blocks.
        v = blockVariance(w(:));
    end
    if v>0, variances(k) = exp(log(v)+2*shift); end
end
var_low = sum(variances);
var_up = sum(sqrt(variances))^2;
end

function v = blockVariance(C)
N = numel(C);
if N<2, v = 0; return; end
sigmaC = var(C,0);
if sigmaC==0, v = 0; return; end
Reff = 1;
for b = 2:max(1,round(N/4))
    nBlock = N-b+1;
    sums = [0;cumsum(C(:))];
    B = (sums(b+1:end)-sums(1:end-b))/b;
    sigmaB = var(B,0)*(nBlock-1)/nBlock;
    Reff = max(Reff,b*sigmaB/sigmaC);
end
v = Reff*sigmaC/N;
end
