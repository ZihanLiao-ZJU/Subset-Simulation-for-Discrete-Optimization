function Z=DiscreteMap(u,xmin,xmax,uniqueValues,mapping)
% Map standard-normal columns to bounded integer designs.
% quantile preserves the original rounded/auction map. uniform gives equal
% integer bins and, for a common all-different range, draws without replacement.
if nargin<4,uniqueValues=true;end
if nargin<5,mapping='auto';end
validateattributes(xmin,{'numeric'},{'real','vector','finite','integer','nonempty'});
validateattributes(xmax,{'numeric'},{'real','vector','finite','integer','nonempty'});
lo=double(xmin(:));hi=double(xmax(:));d=numel(lo);
if numel(hi)~=d || any(lo>hi)
    error('DiscreteMap:Bounds','Supply one nonempty integer interval per coordinate.');
end
if ~((islogical(uniqueValues)||isnumeric(uniqueValues)) && isscalar(uniqueValues) && ...
        isreal(uniqueValues) && isfinite(uniqueValues) && any(uniqueValues==[0,1]))
    error('DiscreteMap:Unique','Unique must be a scalar logical.');
end
if ~isnumeric(u) || ~isreal(u) || ndims(u)>2 || size(u,1)~=d || any(~isfinite(u(:)))
    error('DiscreteMap:Input','Search samples must be finite real dimension-by-sample columns.');
end
if isstring(mapping) && isscalar(mapping),mapping=char(mapping);end
if ~ischar(mapping) || ~isrow(mapping)
    error('DiscreteMap:Mapping','Choose auto, quantile, rank or uniform.');
end
mapping=lower(mapping);
permutationBounds=all(lo==1) && all(hi==d);
if strcmp(mapping,'auto')
    if uniqueValues && permutationBounds,mapping='rank';else,mapping='quantile';end
end
switch mapping
    case 'rank'
        if ~uniqueValues || ~permutationBounds
            error('DiscreteMap:Rank','Rank mapping requires all-different bounds 1:dimension.');
        end
        [~,Z]=sort(u,1);
    case 'quantile'
        U=normcdf(u);U=min(max(U,eps),1-eps);
        X=lo+(hi-lo).*U;
        if uniqueValues
            Z=cont2uniqint(X,lo,hi);
        else
            Z=round(X);Z=min(max(Z,lo),hi);
        end
    case 'uniform'
        U=normcdf(u);U=min(max(U,eps),1-eps);
        if ~uniqueValues
            Z=min(lo+floor((hi-lo+1).*U),hi);
        else
            if ~all(lo==lo(1)) || ~all(hi==hi(1)) || hi(1)-lo(1)+1<d
                error('DiscreteMap:UniformBounds', ...
                    'Uniform all-different mapping requires a common interval containing enough integers.');
            end
            Z=zeros(size(u));
            for k=1:size(u,2)
                pool=lo(1):hi(1);
                for j=1:d
                    index=min(1+floor(U(j,k)*numel(pool)),numel(pool));
                    Z(j,k)=pool(index);pool(index)=[];
                end
            end
        end
    otherwise
        error('DiscreteMap:Mapping','Choose auto, quantile, rank or uniform.');
end
end
