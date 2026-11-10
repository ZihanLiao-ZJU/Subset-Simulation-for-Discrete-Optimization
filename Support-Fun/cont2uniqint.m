function Z=cont2uniqint(X,xmin,xmax)
% Closest all-different integer assignment by the original epsilon auction.
% Feasibility is checked before bidding; impossible intervals never hang.
validateattributes(xmin,{'numeric'},{'real','vector','finite','integer','nonempty'});
validateattributes(xmax,{'numeric'},{'real','vector','finite','integer','nonempty'});
xmin=double(xmin(:));xmax=double(xmax(:));d=numel(xmin);
if numel(xmax)~=d || any(xmin>xmax)
    error('cont2uniqint:Bounds','Supply one nonempty integer interval per coordinate.');
end
if ~isnumeric(X) || ~isreal(X) || ndims(X)>2 || size(X,1)~=d || any(~isfinite(X(:)))
    error('cont2uniqint:Input','Continuous samples must have one finite real row per coordinate.');
end
% Earliest-deadline assignment is an exact feasibility test for intervals.
[~,order]=sort(xmax,'ascend');used=zeros(1,0);
for j=order.'
    candidate=xmin(j);
    while any(used==candidate)
        next=candidate+1;
        if next==candidate,candidate=inf;break;end
        candidate=next;
    end
    if candidate>xmax(j)
        error('cont2uniqint:Infeasible','The intervals admit no all-different integer design.');
    end
    used(end+1)=candidate; %#ok<AGROW>
end
Nsam=size(X,2);Z=zeros(d,Nsam);
if Nsam==0,return;end
S=[];
for i=1:d,S=[S,xmin(i):xmax(i)];end %#ok<AGROW>
S=unique(S,'sorted');M=numel(S);
for k=1:Nsam
    C=abs(X(:,k)-S).^2;
    for i=1:d,C(i,S<xmin(i)|S>xmax(i))=inf;end
    finiteC=C(isfinite(C));base=max(finiteC(:));
    V=base-C;V(~isfinite(C))=-inf;
    epsilon=1e-6;prices=zeros(1,M);r2c=zeros(1,d);c2r=zeros(1,M);U=1:d;
    bids=0;maxBids=max(10000,100*d*M);
    while ~isempty(U)
        bids=bids+1;
        if bids>maxBids
            error('cont2uniqint:AuctionLimit','Integer assignment exceeded its finite bidding limit.');
        end
        i=U(end);U(end)=[];
        [values,index]=sort(V(i,:)-prices,'descend');
        first=find(isfinite(values),1,'first');
        if isempty(first)
            error('cont2uniqint:Assignment','No finite admissible auction candidate remains.');
        end
        j1=index(first);a1=values(first);j2v=-inf;
        if first<numel(values) && isfinite(values(first+1)),j2v=values(first+1);end
        increment=a1-j2v;
        if ~isfinite(increment)||increment<=0,increment=epsilon;else,increment=increment+epsilon;end
        prices(j1)=prices(j1)+increment;
        previous=c2r(j1);c2r(j1)=i;r2c(i)=j1;
        if previous~=0,r2c(previous)=0;U(end+1)=previous;end %#ok<AGROW>
    end
    Z(:,k)=S(r2c).';
    if any(Z(:,k)<xmin|Z(:,k)>xmax) || numel(unique(Z(:,k)))~=d
        error('cont2uniqint:Assignment','Auction returned an invalid all-different design.');
    end
end
end
