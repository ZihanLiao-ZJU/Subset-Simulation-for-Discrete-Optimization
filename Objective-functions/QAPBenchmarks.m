classdef QAPBenchmarks
    % Native permutation benchmarks from QAPLIB, with unmodified data.
    % Cost = sum_ij A(i,j)*B(permutation(i),permutation(j)).
    % https://coral.ise.lehigh.edu/data-sets/qaplib/
    properties (SetAccess=private)
        Name
        Ndim
        Nfun = 1
        Ncst = 0
        Xmin
        Xmax
        Unique = true
        PermutationInvariant = false
        Mapping = 'rank'
        A
        B
    end
    methods
        function obj=QAPBenchmarks(name)
            obj.Name=validatestring(lower(char(name)),{'nug20','nug30','sko49','sko100a'});
            path=fullfile(fileparts(mfilename('fullpath')),'QAPLIB',[obj.Name '.dat']);
            id=fopen(path,'r');
            if id<0,error('QAPBenchmarks:Data','Cannot read instance: %s',path);end
            cleanup=onCleanup(@()fclose(id));
            values=fscanf(id,'%f');n=values(1);
            if n~=fix(n)||n<2||numel(values)~=1+2*n*n||any(~isfinite(values))
                error('QAPBenchmarks:Data','Invalid QAPLIB data in %s.',path);
            end
            obj.Ndim=n;
            obj.A=reshape(values(2:1+n*n),n,n).';
            obj.B=reshape(values(2+n*n:end),n,n).';
            obj.Xmin=ones(1,n);obj.Xmax=n*ones(1,n);
        end

        function F=EvlF(obj,Z)
            if ~isnumeric(Z)||~isreal(Z)||size(Z,1)~=obj.Ndim|| ...
                    any(~isfinite(Z(:)))||any(Z(:)~=fix(Z(:)))|| ...
                    any(Z(:)<1)||any(Z(:)>obj.Ndim)|| ...
                    any(any(diff(sort(Z,1),1,1)==0))
                error('QAPBenchmarks:Permutation','Expected permutation columns of 1:Ndim.');
            end
            Z=double(Z);F=zeros(1,size(Z,2));
            for k=1:size(Z,2)
                p=Z(:,k);F(k)=sum(obj.A.*obj.B(p,p),'all');
            end
        end

        function C=EvlCST(~,Z)
            C=zeros(0,size(Z,2));
        end

        function reference=Oracle(obj)
            % Reference metadata is for external validation only. Neither the
            % optimizer nor the driver reads it; no solution is supplied.
            switch obj.Name
                case 'nug20',value=2570;kind='published proven optimum';
                case 'nug30',value=6124;kind='published proven optimum';
                case 'sko49',value=23386;kind='QAPLIB listed best known feasible value';
                case 'sko100a',value=152002;kind='QAPLIB listed best known feasible value';
            end
            reference=struct('name',obj.Name,'kind',kind,'objectives',value, ...
                'candidate_count',factorial(obj.Ndim),'representatives',[], ...
                'source','https://coral.ise.lehigh.edu/data-sets/qaplib/qaplib-problem-instances-and-solutions/');
        end
    end
end
