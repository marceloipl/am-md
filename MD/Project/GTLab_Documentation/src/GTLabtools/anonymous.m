disp("List of useful anonymous functions");

nlin=@(x) size(x,1);
ncol=@(x) size(x,2);
idn=@(x) (1:x)';
id=@(x) (1:numel(x))';

% collection of divisors of a given number
div=@(n) find(mod(n,1:n)==0);
% and respective Hasse diagram
u=div(1000); [x,y]=ndgrid(id(u));
R=(mod(u(y),u(x))==0); 


% properties of relations
isreflexive=@(R) all(diag(R));
istransitive=@(R) all(all(logical(R*R)<=R));
issymmetric=@(R) isequal(R,R');
isantisymmetric=@(R) all(all((R&R')<=eye(size(R))));


% the composition of endomaps as if they were row vectors
composerows=@(A,J) A(sub2ind(size(J),repmat((1:size(J,1))',1,size(J,2)),J)); % if AJ=composerows(A,J) then, for each i, AJ(i,:)=A(i,J(i,:))
% the bijection between A and B^X
s2i=@(nB,s) (s-1)*(nB.^((1:size(s,2))'-1))+1; % nX=size(s,2)
i2s=@(i,nB,nX) mod(floor((i(:)-1)./nB.^((1:nX)-1)),nB)+1; % using the fact that i./s expands i to the size of si
circ=@(i,j,nB,nX,nZ) s2i(composerows(i2s(i,nB,nX),i2s(j,nX,nZ)),nB); % compose (i,nB,nX) with (j,nX,nZ), where (i,nB,nX) represents the morphism X-->B indexed by i in the ser B^X


% random partition (Aluno Paulo Rodrigues --- 2223323)
partition=@(n,k) diff([0 sort(randperm(n-1,k-1)) n]); 
% returns a vector p with k random integers ranging in 1:(n-(k-1)) such that sum(p)=n

% signature interpolation
siginterp=@(sig,n) diff([0, round(cumsum(sig)/sum(sig)*n)]);

% signature to connected components and back
sig2conn=@(sig) repelem(1:numel(sig),sig)
conn2sig=@(conn) full(sparse(conn,1,1))' % histc(conn,unique(conn))

% signature to cycles and back
sig2cycles=@(sig) 1+idn(sum(sig))+sparse(cumsum(sig),1,-sig);
cycles2sig=@(phi) conn2sig(orbits(phi))



%%%----------- To be sorted out
%{
%% The functor from Finord to Mat
fin2mat=@(n,v) sparse(id(v),v,1,numel(v),n);

%mat2fin=@(A) [k,v]=max(logical(A),[],2);

function v=mat2fin(A)
[k,p]=max(logical(A'));
v=p(k)';
end



% given p returns the shift so that
% phi=2:sum(p)+1;
% phi(cumsum(p))=shift(p);

shift=@(p) sort(mod(cumsum(sort(p))+1,sum(p)));
sigshift=@(p) sort(mod(cumsum(sort(p))+1,sum(p)));
% not to colide with Matlab shift.m

% p=sort(partition(12,3))
p=[4 5 6]
phi=2:sum(p)+1;
phi(cumsum(p))=shift(p);

% or
sig2phi=@(sig) cumsum((sparse([cumsum(sig) shift(sig)],1,[-sig shift(sig)-1],sum(sig),1))+1)+1;

% and retuns connected components
sig2conn=@(p) full(cumsum(sparse(shift(sort(p)),1,1,sum(p),1)))';

% conversely
conn2sig=@(conn) diff([0 find(diff(conn)) numel(conn)]);

% or
sig2conn=@(sig) repelem(1:numel(sig),sig)
conn2sig=@(conn) histc(conn,unique(conn))


% to interpolate
nf=10
f=sig2conn([4 5 6])

r=floor(find(diff(f))/numel(f)*nf)
fnew=zeros(1,nf)
fnew([1 r+1])=1
f_=cumsum(fnew)

interpolate=@(conn,m) full(cumsum(sparse([1 ...
floor(find(diff(conn))/numel(conn)*m)...
+1],1,1,m,1)))'-1;

%}