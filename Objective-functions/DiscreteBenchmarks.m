classdef DiscreteBenchmarks
    % Classical discrete families and integer-coded grids of test functions.
    % All objectives are minimized. The knapsack and geometric instances are
    % constructed here; they are not named instances from a benchmark archive.
    % LOTZ/COCZ definitions: Laumanns, Thiele and Zitzler (2004),
    % https://sop.tik.ee.ethz.ch/publicationListFiles/ltz2004a.pdf
    % Grid objectives retain the continuous formulas but restrict their domain;
    % they are not claimed to be original, published discrete instances.
    % ZDT1: https://sop.tik.ee.ethz.ch/publicationListFiles/zdt2000a.pdf
    % DTLZ2: https://sop.tik.ee.ethz.ch/publicationListFiles/dtlz2002a.pdf
    properties (SetAccess=private)
        Name
        Ndim
        Nfun = 1
        Ncst = 0
        Xmin
        Xmax
        Unique = false
        PermutationInvariant = false
        Mapping = 'quantile'
        Levels = []
        GridLower = []
        GridUpper = []
        Weights = [2,3,6,7,5,9,4,8,1,10]
        Profits = [6,5,8,9,6,12,7,10,2,13]
        Capacity = 25
        Coordinates = []
        Distances = []
    end
    methods
        function obj=DiscreteBenchmarks(name,varargin)
            if nargin<1,name='onemax';end
            p=inputParser;
            p.addParameter('Dimension',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=1&&x==fix(x)));
            p.addParameter('Levels',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=2&&x==fix(x)));
            p.parse(varargin{:});opt=p.Results;
            obj.Name=validatestring(lower(char(name)), ...
                {'onemax','knapsack','pmedian','tsp','lotz','cocz','biknapsack', ...
                'grid_rastrigin','grid_rosenbrock','grid_ackley','grid_zdt1','grid_dtlz2'});
            grid=strncmp(obj.Name,'grid_',5);
            if ~grid && ~isempty(opt.Levels)
                error('DiscreteBenchmarks:Levels','Levels applies only to integer-coded grid examples.');
            end
            if ~isempty(opt.Dimension) && ~grid && ...
                    ~any(strcmp(obj.Name,{'onemax','lotz','cocz'}))
                error('DiscreteBenchmarks:Dimension','This fixed instance has a prescribed dimension.');
            end
            switch obj.Name
                case 'onemax'
                    obj.Ndim=20;
                case {'knapsack','biknapsack'}
                    obj.Ndim=numel(obj.Weights);
                    if strcmp(obj.Name,'biknapsack'),obj.Nfun=2;obj.Ncst=1;end
                case 'pmedian'
                    obj.Ndim=3;obj.Unique=true;obj.PermutationInvariant=true;
                    [xx,yy]=ndgrid(0:3,0:2);obj.Coordinates=[xx(:).';yy(:).'];
                    obj.Distances=DiscreteBenchmarks.distances(obj.Coordinates);
                case 'tsp'
                    obj.Ndim=8;obj.Unique=true;obj.Mapping='rank';
                    theta=2*pi*(0:7)/8;obj.Coordinates=[cos(theta);sin(theta)];
                    obj.Distances=DiscreteBenchmarks.distances(obj.Coordinates);
                case {'lotz','cocz'}
                    obj.Ndim=12;obj.Nfun=2;
                case {'grid_rastrigin','grid_rosenbrock','grid_ackley'}
                    obj.Ndim=20;
                case 'grid_zdt1'
                    obj.Ndim=30;obj.Nfun=2;
                case 'grid_dtlz2'
                    obj.Ndim=12;obj.Nfun=3;
            end
            if ~isempty(opt.Dimension),obj.Ndim=opt.Dimension;end
            if strcmp(obj.Name,'cocz') && mod(obj.Ndim,2)
                error('DiscreteBenchmarks:Dimension','COCZ requires an even dimension.');
            end
            if any(strcmp(obj.Name,{'grid_rosenbrock','grid_zdt1'})) && obj.Ndim<2
                error('DiscreteBenchmarks:Dimension','This objective requires at least two variables.');
            elseif strcmp(obj.Name,'grid_dtlz2') && obj.Ndim<3
                error('DiscreteBenchmarks:Dimension','Three-objective DTLZ2 requires at least three variables.');
            end
            obj.Xmin=zeros(1,obj.Ndim);obj.Xmax=ones(1,obj.Ndim);
            if grid
                obj.Mapping='uniform';
                obj.Levels=101;
                if ~isempty(opt.Levels),obj.Levels=opt.Levels;end
                if mod(obj.Levels,2)==0
                    error('DiscreteBenchmarks:Levels','Grid examples require odd Levels so their known minimizer is represented.');
                end
                if strcmp(obj.Name,'grid_rosenbrock') && mod(obj.Levels-1,4)
                    error('DiscreteBenchmarks:Levels','Rosenbrock on [-2,2] requires Levels-1 divisible by four to represent x=1.');
                end
                switch obj.Name
                    case 'grid_rastrigin',obj.GridLower=-5;obj.GridUpper=5;
                    case 'grid_rosenbrock',obj.GridLower=-2;obj.GridUpper=2;
                    case 'grid_ackley',obj.GridLower=-32;obj.GridUpper=32;
                    otherwise,obj.GridLower=0;obj.GridUpper=1;
                end
                obj.Xmax(:)=obj.Levels-1;
            elseif strcmp(obj.Name,'pmedian')
                obj.Xmin(:)=1;obj.Xmax(:)=size(obj.Coordinates,2);
            elseif strcmp(obj.Name,'tsp')
                obj.Xmin(:)=1;obj.Xmax(:)=obj.Ndim;
            end
        end

        function F=EvlF(obj,Z)
            % One design per column; no random numbers or hidden model state.
            if size(Z,1)~=obj.Ndim,error('DiscreteBenchmarks:Shape','Expected Ndim rows.');end
            switch obj.Name
                case 'onemax'
                    F=obj.Ndim-sum(Z,1);
                case 'knapsack'
                    % Integer excess >=1: this penalty makes every infeasible
                    % value positive, while every feasible value is <=0.
                    excess=max(0,obj.Weights*Z-obj.Capacity);
                    F=-obj.Profits*Z+(sum(obj.Profits)+1)*excess;
                case 'pmedian'
                    d=obj.Distances(:,Z(:));
                    d=reshape(d,size(obj.Distances,1),obj.Ndim,size(Z,2));
                    F=reshape(sum(min(d,[],2),1),1,[]);
                case 'tsp'
                    following=Z([2:end,1],:);
                    edges=sub2ind(size(obj.Distances),Z,following);
                    F=sum(obj.Distances(edges),1);
                case 'lotz'
                    F=-[sum(cumprod(Z,1),1);sum(cumprod(flipud(1-Z),1),1)];
                case 'cocz'
                    h=obj.Ndim/2;
                    F=-[sum(Z,1);sum(Z(1:h,:),1)+sum(1-Z(h+1:end,:),1)];
                case 'biknapsack'
                    F=[-obj.Profits*Z;obj.Weights*Z];
                case 'grid_rastrigin'
                    X=obj.gridValues(Z);F=sum(X.^2+10*(1-cos(2*pi*X)),1);
                case 'grid_rosenbrock'
                    X=obj.gridValues(Z);
                    F=sum(100*(X(2:end,:)-X(1:end-1,:).^2).^2+(X(1:end-1,:)-1).^2,1);
                case 'grid_ackley'
                    X=obj.gridValues(Z);
                    % Equivalent to the standard Ackley expression, with
                    % expm1 preserving an exactly zero value at the origin.
                    F=-20*expm1(-.2*sqrt(mean(X.^2,1))) ...
                        -exp(1)*expm1(mean(cos(2*pi*X),1)-1);
                case 'grid_zdt1'
                    X=obj.gridValues(Z);g=1+9*mean(X(2:end,:),1);
                    F=[X(1,:);g.*(1-sqrt(X(1,:)./g))];
                case 'grid_dtlz2'
                    X=obj.gridValues(Z);g=sum((X(3:end,:)-.5).^2,1);
                    c=cos(pi*X(1:2,:)/2);s=sin(pi*X(1:2,:)/2);
                    c(X(1:2,:)==1)=0;s(X(1:2,:)==0)=0;
                    F=(1+g).*[c(1,:).*c(2,:);c(1,:).*s(2,:);s(1,:)];
            end
        end

        function C=EvlCST(obj,Z)
            if obj.Ncst,C=obj.Weights*Z-obj.Capacity;
            else,C=zeros(0,size(Z,2));end
        end

        function reference=Oracle(obj)
            % Reference calculation for external validation, never called by
            % drive_discrete or by an optimization model/sampler.
            switch obj.Name
                case 'onemax'
                    Z=ones(obj.Ndim,1);F=0;count=2^obj.Ndim;kind='analytic';
                case 'lotz'
                    k=0:obj.Ndim;Z=double((1:obj.Ndim).'<=k);
                    F=-[k;obj.Ndim-k];count=2^obj.Ndim;kind='analytic Pareto front';
                case 'cocz'
                    h=obj.Ndim/2;k=0:h;
                    Z=[ones(h,h+1);double((1:h).'<=k)];
                    F=-[h+k;obj.Ndim-k];count=2^obj.Ndim;kind='analytic Pareto front';
                case {'knapsack','biknapsack'}
                    codes=uint16(0:2^obj.Ndim-1);Z=zeros(obj.Ndim,numel(codes));
                    for bit=1:obj.Ndim,Z(bit,:)=bitget(codes,bit);end
                    weight=obj.Weights*Z;profit=obj.Profits*Z;
                    keep=weight<=obj.Capacity;Z=Z(:,keep);
                    if strcmp(obj.Name,'knapsack')
                        [value,index]=max(profit(keep));Z=Z(:,index);F=-value;
                    else
                        F=[-profit(keep);weight(keep)];
                        [F,index]=unique(F.','rows','stable');F=F.';Z=Z(:,index);
                        keep=DiscreteBenchmarks.pareto(F);Z=Z(:,keep);F=F(:,keep);
                    end
                    count=2^obj.Ndim;kind='exhaustive binary enumeration';
                case 'pmedian'
                    Z=nchoosek(1:size(obj.Coordinates,2),obj.Ndim).';count=size(Z,2);
                    % Independent direct loop over demands and candidate sets.
                    values=zeros(1,count);
                    for j=1:count
                        for demand=1:size(obj.Coordinates,2)
                            delta=obj.Coordinates(:,Z(:,j))-obj.Coordinates(:,demand);
                            values(j)=values(j)+min(sqrt(sum(delta.^2,1)));
                        end
                    end
                    [F,index]=min(values);Z=Z(:,index);kind='exhaustive subset enumeration';
                case 'tsp'
                    Z=[ones(1,factorial(obj.Ndim-1));perms(2:obj.Ndim).'];
                    count=size(Z,2);values=zeros(1,count);
                    for j=1:count
                        tour=obj.Coordinates(:,Z([1:end,1],j));
                        values(j)=sum(sqrt(sum(diff(tour,1,2).^2,1)));
                    end
                    [F,index]=min(values);Z=Z(:,index);kind='exhaustive rooted-tour enumeration';
                case {'grid_rastrigin','grid_rosenbrock','grid_ackley'}
                    code=(obj.Levels-1)/2;
                    if strcmp(obj.Name,'grid_rosenbrock'),code=3*(obj.Levels-1)/4;end
                    Z=repmat(code,obj.Ndim,1);F=0;
                    count=double(obj.Levels)^obj.Ndim;kind='analytic grid minimum';
                case 'grid_zdt1'
                    t=(0:obj.Levels-1)/(obj.Levels-1);
                    Z=zeros(obj.Ndim,obj.Levels);Z(1,:)=0:obj.Levels-1;
                    F=[t;1-sqrt(t)];count=double(obj.Levels)^obj.Ndim;
                    kind='analytic complete grid Pareto front';
                case 'grid_dtlz2'
                    % At g=0, all angular grid points are nondominated because
                    % their nonnegative objective vectors have unit norm.
                    % The north pole is independent of the second angle;
                    % retain one representative to avoid duplicate PF points.
                    [a,b]=ndgrid(0:obj.Levels-2,0:obj.Levels-1);
                    angles=[a(:).';b(:).'];angles(:,end+1)=[obj.Levels-1;0];
                    Z=repmat((obj.Levels-1)/2,obj.Ndim,size(angles,2));Z(1:2,:)=angles;
                    t=angles/(obj.Levels-1);c=cos(pi*t/2);s=sin(pi*t/2);
                    c(t==1)=0;s(t==0)=0;
                    F=[c(1,:).*c(2,:);c(1,:).*s(2,:);s(1,:)];
                    count=double(obj.Levels)^obj.Ndim;kind='analytic complete grid Pareto front';
            end
            reference=struct('solutions',Z,'objectives',F,'candidate_count',count, ...
                'kind',kind,'tolerance',1e-10);
        end
    end
    methods (Access=private)
        function X=gridValues(obj,Z)
            % The optimizer still sees integer codes, never real-valued genes.
            X=obj.GridLower+(obj.GridUpper-obj.GridLower)*Z/(obj.Levels-1);
        end
    end
    methods (Static,Access=private)
        function D=distances(x)
            D=sqrt(sum((reshape(x,size(x,1),[],1)-reshape(x,size(x,1),1,[])).^2,1));
            D=reshape(D,size(x,2),size(x,2));
        end
        function mask=pareto(F)
            % Deliberately independent of the algorithm's NonDomRnk helper.
            mask=true(1,size(F,2));
            for j=1:size(F,2)
                mask(j)=~any(all(F<=F(:,j),1)&any(F<F(:,j),1));
            end
        end
    end
end
