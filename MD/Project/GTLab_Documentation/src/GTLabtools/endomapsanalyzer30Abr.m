function [Rstar,p,q]=endomapsanalyzer30Abr(k)
%  [Rstar,p,q]=endomapsanalyzer(k)
% New version 5Jun2025 -- runs with new functions hasse.m and analyzerelation.m, no longer using the auxiliar functions on this file
% 30Abr2025; 13:00 - with small bugs fixed in order to run in Matlab
% 
% Possible solution to Class Activity proposed in Guide4Classes.pdf version 3

version='05Jun2025'

switch version
case '05Jun2025'
% The bijection between A and B^X with parameters
s2i=@(nB,s) (s-1)*(nB.^((1:size(s,2))'-1))+1; % nX=size(s,2)
i2s=@(i,nB,nX) mod(floor((i(:)-1)./nB.^((1:nX)-1)),nB)+1; % using the fact that i./s expands i to the size of s

nB=3, nX=3, nA=nB^nX

% The bijection between A and B^X for the specific case nX, nB and nA
s2i=@(s) s2i(nB,s); 
i2s=@(i) i2s(i,nB,nX);
% checking bijection between A and B^X
isequal((1:nA)',s2i(i2s(1:nA)))
isequal(i2s(1:nA),i2s(s2i(i2s(1:nA))))


% the composition of endomaps as if they were row vectors
composerows=@(A,J) A(sub2ind(size(J),repmat((1:size(J,1))',1,size(J,2)),J)); % if AJ=composerows(A,J) then, for each i, AJ(i,:)=A(i,J(i,:))

Linv=all(sort(i2s(1:nA),2)==repmat(1:nX,nA,1),2); % invertible maps
Lsort=all(diff(i2s(1:nA)')'>=0,2); % monotone maps


% Each one of the six cases
switch k
case 1 % Linv, d inclusion, c inverse
L=Linv;
d=find(L);
% to compute c we observe:
for i=1:nnz(L)
invi=zeros(1,nX);
invi(i2s(d(i)))=(1:nX);
c(i)=s2i(invi);
end
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(d,1,d,nA,1)) full(sparse(d,1,c,nA,1))]

case 2 % Linv, d trivial, c inclusion
L=Linv;
% the trivial endomap is i2s(22) or [1 2 3]:
d=repmat(s2i(1:nX),1,nnz(L))
c=find(L);
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(c,1,d,nA,1)) full(sparse(c,1,c,nA,1))]

case 3 % Linv, d inclusion, c trivial
L=Linv;
d=find(L);
c=repmat(s2i(1:nX),1,nnz(L))
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(d,1,d,nA,1)) full(sparse(d,1,c,nA,1))]

case 4 % Lsort, d trivial, c inclusion
L=Lsort;
c=find(L);
d=repmat(s2i(1:nX),1,nnz(L))
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(c,1,d,nA,1)) full(sparse(c,1,c,nA,1))]

case 5 %Lsort, d inclusion, c trivial
L=Lsort;
d=find(L);
c=repmat(s2i(1:nX),1,nnz(L))
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(d,1,d,nA,1)) full(sparse(d,1,c,nA,1))]

case 6 %Lsort*Lsort, d first projection, c second projection
L=Lsort;

[v1,v2]=ndgrid(find(L));
d=v1(:);
c=v2(:);
% check for accurary
[(1:nnz(L)^2)', d c]

end % switch

% defining and analyzing the Simplicial Relation

% given L logical vector of size nA, together with d and c
% the simplicial relation is
[a i]=ndgrid(1:nA,1:numel(d));
b=s2i(composerows(i2s(d(i)), composerows(i2s(a),i2s(c(i)))));
S=logical(sparse(a,b,1,nA,nA));
%[a,b]=find(S);

[Rstar,p,q]=analyzerelation(S);
disp('Analyzing the smallest equivalence relation containing S, encoded as an endomap where each element points to its representative:')
q(:)'
linkdraw(homology(q),q);
pause
SM=sortrows([q, (1:nA)' i2s(1:nA)]);
showmethemaps(SM);
pause
subplot(1,1,1)
disp('Analyzing the partial order induced by S, showing the respective Hasse diagram and adjacency matrix')
u=unique(p);
Rfull=full(Rstar(u,u))
[x,y]=hasse(Rfull,u);


case '30Abr2025'  % old version
% setting the stage
nX=3, nB=nX, nA=nB^nX
i2s=@(i) mod(floor(bsxfun(@rdivide,i(:)-1,nB.^((1:nX)-1))),nB)+1;
s2i=@(s) (s-1)*(nB.^((1:nX)'-1))+1;
%A=(1:nA)';
% checking bijection between A and B^X
isequal((1:nA)',s2i(i2s(1:nA)))
isequal(i2s(1:nA),i2s(s2i(i2s(1:nA))))

% the composition of endomaps as if they were row vectors
composerows=@(A,J) A(sub2ind(size(J),repmat((1:size(J,1))',1,size(J,2)),J)); % if AJ=composerows(A,J) then, for each i, AJ(i,:)=A(i,J(i,:))

Linv=all(sort(i2s(1:nA),2)==repmat(1:nX,nA,1),2); % invertible maps
Lsort=all(diff(i2s(1:nA)')'>=0,2); % monotone maps


% Each one of the six cases
switch k
case 1 % Linv, d inclusion, c inverse
L=Linv;
d=find(L);
% to compute c we observe:
for i=1:nnz(L)
invi=zeros(1,nX);
invi(i2s(d(i)))=(1:nX);
c(i)=s2i(invi);
end
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(d,1,d,nA,1)) full(sparse(d,1,c,nA,1))]

case 2 % Linv, d trivial, c inclusion
L=Linv;
% the trivial endomap is i2s(22) or [1 2 3]:
d=repmat(s2i(1:nX),1,nnz(L))
c=find(L);
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(c,1,d,nA,1)) full(sparse(c,1,c,nA,1))]

case 3 % Linv, d inclusion, c trivial
L=Linv;
d=find(L);
c=repmat(s2i(1:nX),1,nnz(L))
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(d,1,d,nA,1)) full(sparse(d,1,c,nA,1))]

case 4 % Lsort, d trivial, c inclusion
L=Lsort;
c=find(L);
d=repmat(s2i(1:nX),1,nnz(L))
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(c,1,d,nA,1)) full(sparse(c,1,c,nA,1))]

case 5 %Lsort, d inclusion, c trivial
L=Lsort;
d=find(L);
c=repmat(s2i(1:nX),1,nnz(L))
% check for accurary
[(1:nA)', i2s(1:nA) L full(sparse(d,1,d,nA,1)) full(sparse(d,1,c,nA,1))]

case 6 %Lsort*Lsort, d first projection, c second projection
L=Lsort;

[v1,v2]=ndgrid(find(L));
d=v1(:);
c=v2(:);
% check for accurary
[(1:nnz(L)^2)', d c]

end % switch

% defining and analyzing the Simplicial Relation

% given L logical vector of size nA, together with d and c
% the simplicial relation is
[a i]=ndgrid(1:nA,1:numel(d));
b=s2i(composerows(i2s(d(i)), composerows(i2s(a),i2s(c(i)))));
S=logical(sparse(a,b,1,nA,nA));
%[a,b]=find(S);

[Rstar,p,q]=analyzerelation(S);
% analyzing the smallest equivalence relation containing S
linkdraw(homology(q),q);
pause
SM=sortrows([q, (1:nA)' i2s(1:nA)]);
showmethemaps(SM);
pause
subplot(1,1,1)
% analyzing the partial order induced by S
u=unique(p)
full(Rstar(u,u))
%analyzerelation(Rstar(u,u));
[x,y]=hasse(Rstar(u,u),u);

end % switch version

end % endomaps main function

% From here on, the auxiliar functions, are no longer being used in version 05Jun2025, for that reason it is now commented as a block
%{
% ----------- Auxiliar functions------------------
function [Rstar,p,q]=analyzerelation(R)
% [Rstar,p,q]=analyzerelation(R)
% 30Abr2025
% Given the adjacency matrix R of a binary relation on the set A, returns Rstar, the reflexive and transitive closure of R, as well as:
% (a) p:A-->P, the partition induced by the equivalence relation Rstar&Rstar', so that the relation Rstar(unique(p),unique(p)) is a partial order
% (b) q:A-->Q, the partition induced by the smallest equivalence relation, E, containing R, so that, E can be obtained as:
% [i,j]=ndgrid(1:size(R,1),1:size(R,2)); E=(q(i)==q(j));
% 
% In addition it displays the properties reflexive, transitive, etc

sprintf('is reflexive? %d \n',full(all(diag(R))))
sprintf('is transitive? %d \n',full(all(all(logical(R*R)<=R))))
sprintf('is symmetric? %d \n',isequal(R,R'))
sprintf('is anti-symmetric? %d \n',full(all(all((R&R')<=eye(size(R))))))

Rtran=transitiveclosure(R);
Rstar=Rtran | eye(size(R)); % reflexive and transitive closure
[d,c]=find(Rstar);
q=coeq(d,c);
[d,c]=find(Rstar&Rstar');
p=coeq(d,c);


end % analyzerelation
%-----------------------------------------------------
function Rtran=transitiveclosure(R1)
% R=transitiveclosure(R1)
% 30Abr2025
% Rtran is the transitive closure of R1
%
% properties of relations
% isreflexive=@(R) all(diag(R));
% istransitive=@(R) all(all(logical(R*R)<=R));
% issymmetric=@(R) isequal(R,R');
% isantisymmetric=@(R) all(all(R&R'<=eye(size(R))));

R=R1; % Get the transitive closure of R1
Rtran=logical(R+R*R1);
while ~isequal(R,Rtran)
    R=Rtran;
    Rtran=logical(R+R*R1);
end

end % transitiveclosure
%--------------------------------------------------------
function [x,y]=hasse(R,D)
% [x,y]=hasse(R,D)
% 30Abr2025
% given the adjacency matrix R of a partial order with elements D, return the x and y coordinates to draw the respective hasse diagram

nD=numel(D);
if ~isequal(size(R),[nD nD]), error('R must be nD-by-nD, with nD=numel(D)'), end
if ~all(all(logical(R*R)<=R)) | ~all(diag(R)) | ~all(all((R&R')<=eye(size(R))))
    error('R is not a partial order')
end

s=sum(R)'-sum(R,2);
[~,~,q]=unique(s);
[qs,sq]=sort(q);
%disp([s(sq) qs   mod((1:nD)'-1,qs)+1 (1:nD)'])
y(sq)=s(sq);
v=hist(qs,1:max(qs))';
x(sq)=(mod((1:nD)'-1,v(qs))+1) - v(qs)/2  + 0.5*sq/nD;

y=(y-min(y))/(max(y)-min(y));
y(isnan(y))=0;
x=(x-min(x))/(max(x)-min(x));

% in order to draw the lines in the hasse diagram
R2=(R^2==2);
[i,j]=find(R2|eye(size(R)));
du=x(j)-x(i);
dv=y(j)-y(i);
%disp('input to quiver'), disp([x(i),y(i),du,dv])
hq=quiver(x(i),y(i),du,dv,0);
 set(hq,'showarrowhead','off')
text(x,y,num2str(D(:))), axis([min(x) max(x) min(y)-eps max(y)+eps]), axis off
end % hasse
%-----------------------------------------------
function [q,p1,p2]=coeq(d,c)
% [q,p1,p2]=coeq(d,c)
%
% given a graph (d,c), considered as two maps d,c:A-->B with A=(1:nA)' and B=(1:nB)', and such that d is surjective, returns q=coeq(d,c) the categorical coequalizer of d and c, so that q(d)==q(c) and for every  f:B-->C with f(d)==f(c) there exists a unique fbar:(1:numel(u)')-->C such that f=fbar(r(q)) for a unique fbar obtained as fbar=f(r) with [u,r,q]=unique(q); this replaces the output q indexed by representatives into a new q with linear indexes, so that the old q is recovered as u(q).
% the output values p1,p2 are such that if R=sparse(p1,,p2) then R+R' is the kernel pair relation of q 

nB=max(d); B=(1:nB)';
if ~isequal(unique(d),B)
    fprintf('Warning: d is not surjective, q is not guaranteed to be a coequalizer \n');
end

R0=sparse(d,c,1,nB,nB);
R1=triu(R0+R0');
R=R1; % Get the transitive closure of R1
Rtran=logical(R+R*R1);
while ~isequal(R,Rtran)
    R=Rtran;
    Rtran=logical(R+R*R1);
end

[s,t]=find(R+eye(nB));
E=[s,t]; Erows=sortrows(E);
p1=Erows(:,1); p2=Erows(:,2);
q(p1)=p2; q=q(:);

end % of function coeq
%---------------------------------------
function [g,Hc,Hnc]=homology(phi)
% [g,Hc,Hnc]=homology(phi)
%
% 24May2024 - simplified and tested (the old version is homology_.m)
% given an endomap phi:X-->X with X=(1:numel(phi)) computes the orbits of phi and returns a complex vector g of the same size as phi with geometic realization together with Hc an homology characteristic vector with the number of elements is each cyclic component as well as Hnc a homology matrix with information on each non-cyclic connected component
%
% is an alternative to orbits24.m, orbits24_inprogress and to homologymatrix.m (and to homology_.m)
%
% Example (see also at the end of file)
% phi=[2 5 5, 15 15 16, 9 9 10, 11 8 13, 14 12 19, 21 16 19, 19 21 17, 17 ]
% [g,Hc,Hnc]=homology(phi)
% linkdraw23(phi,g,1:numel(phi))
%
% history: last stabe version24May; bug fixed on line 33 (31May2024)

[orb,ord,psi,deg,init,term,prin,conn]=orbits(phi);
idA=(1:numel(phi))';
g=zeros(size(phi)); % initialize g, the geometrical realization of phi

% the case of cyclic components
iscyc=(deg==-1);
cyc1=find(iscyc & ord==1); % indexes for first cycles
Hc=sortrows([ord(psi(cyc1)), orb(cyc1)]); % the second row contains the reference to the orbit, should be discarded for homology comparisons
% Example: >> phi=[2 3 4 1 6 7 5 9 8]; 
% Hc =[ 2   3   4]'

% the geometric realization for the cyclic components
SM=sortrows([orb(iscyc), ord(iscyc), idA(iscyc)]);
g(SM(:,3))=(1i*2*pi*linspace(1/nnz(iscyc),1,nnz(iscyc)));


% the case of non-cyclic components
connop=conn'; prinop=prin'; % fixing a bug 31May2024
SM=sortrows([connop(orb(~iscyc)), prinop(orb(~iscyc)), ord(phi(term(orb(~iscyc)))), orb(~iscyc), ord(~iscyc), idA(~iscyc)]);  % sorting the data to be used further on
[~,r,q]=unique(conn); 
% conn is replaced by q with linear indexes and r(q)=conn'
Prin=sparse(q,prin,ord(term)); % combines each index of connected components with each one of its principal secondary tertiary components etc

Hnc=zeros(numel(r),6); Hnc(:,6)=r; % the 6th column is only for reference
% the total number of points in each unique connected component
Hnc(:,1)=full(sum(Prin,2));
% the period of each unique connected component
Hnc(:,2)=ord(term(r))-ord(phi(term(r)))+1;
% the degrees of each connected component
for i=1:numel(r)
isconni=(conn(orb(~iscyc))==r(i))';
Hnc(i,3)=nnz(isconni & deg(~iscyc)==0); % initial points
Hnc(i,4)=nnz(isconni & deg(~iscyc)==2); % junction points with deg==2
Hnc(i,5)=nnz(isconni & deg(~iscyc)>=3); %junction points with deg>=3
end

% computing the starting and ending angles for each unique connected component
th=2*pi*(cumsum([0; Hnc(:,1)]))/sum(Hnc(:,1)); 
% each component starts at th(q(i)) and ends and th(q(i)+1) angles between 0 and 2pi

% computing the geometry of each unique connected component associated with its principal, secondary etc, sub-components
[I,J]=find(Prin); 
% each pair (I(i),J(i)) is an index of a connected component (I(i)) combined with an index of a principal, secondary etc sub-component J(i); the component indexed by I(i) is r(I(i)) 
for i=1:numel(I)
    iscompi=(SM(:,1)==r(I(i)) & SM(:,2)==J(i));
        ni=nnz(iscompi);
        thi=th(I(i));
        rangi= 0.5 *(1 + linspace(-1+1/ni,1-1/ni,ni));
        rangethi=min(pi/2,th(I(i)+1)-thi);
        g(SM(iscompi,6))=log(sqrt(J(i)+1)) + 1i*(thi+ rangi*rangethi);
end

g=exp(g);

Hnc=sortrows(Hnc(:,1:5)); % sorting and discarding the 6th column
Hc=Hc(:,1); % removing the 2nd column in Hc which was just a reference to the cyclic component 

end % function homology
%----------------------------------------------------
function [orb, ord, psi, deg, init, term, prin, conn] = orbits(phi)
% [idn, phi, orb, ord, psi, deg, ini, term, prin, conn] = orbits(phi)
%%
% Date: 2Fev2024
% given a vector phi
try
    phi=phi(:); % force phi to be a column vector
    phi(phi); % make sure that phi can be interpreted as an endomap
catch ME
    disp('The input vector phi must be interpreted as an endomap from the set of indexes 1:numel(phi) to itself:')
    error(ME.message)
end
n=numel(phi);
idn=(1:n)'; % the identity vector of the same size as phi
% initializing outputs
orb=zeros(n,1);
ord=zeros(n,1);
psi(phi)=idn; psi=psi(:); % ensure column vector
deg=repmat(-1,n,1);  % the degree is initiated by default as cyclic
% determining initial points (those with degree 0)
init=idn(~ismember(idn,phi))'; % vector with indexes of initial points
m=numel(init); % number of initial elements
term=zeros(1,m); % the vector term will contain the terminal point of each initial point with respect to the choice of the principal components
prin=zeros(1,m); % the vector prin will contain the principal components, secondary, terciary, etc
conn=zeros(1,m); % the connected component, relative to the principal component

% the iterative process
for j=1:numel(init)
    a=init(j); k=1;
    orb(a)=j;
    ord(a)=k;
    psi(a)=a; % the initial points become fixed points when inverted
    deg(a)=0;
    while ~orb(phi(a))   % phi(a) is not yet a junction point
        psi(phi(a))=a; % is defined before a is uptaded as a=phi(a)
        a=phi(a); k=k+1;  % a is updated to phi(a) and k to k+1
        orb(a)=j;
        ord(a)=k;
        deg(a)=1;
    end
    % phi(a) is a junction point (orb(phi(a)) is not zero)
    if orb(a)==orb(phi(a))  %  principal component
        deg(phi(a))=deg(phi(a))+1;
        term(j)=a;
        prin(j)=1;
        conn(j)=j;  % the same as orb(a)
    else % secondary component
        deg(phi(a))=deg(phi(a))+1;
        term(j)=a; % each secondary component will use the terminal point of its principal component
        prin(j)=prin(orb(phi(a)))+1;
        conn(j)=conn(orb(phi(a)));
    end
end

% It remains to compute the cyclic orbits
% e.g. phi=[2 3 4 5 1 1 6 6 1 3 12 13 14 15 12 13 16 17]'; phi=[phi' 20 21 22 19 23 25 24]; disp([idphi phi orb ord psi])
q=m;  % m=numel(init), initialize component label, in oder to continue the non-cyclic components
while any(~orb)
    a=find(~orb,1);
    q=q+1;  % update component lablel
    orb(a)=q;  % set component label
    ord(a)=1;  % initialize counting the order of current component
    %deg(a)=-1;  % we have choosen to single out the cyclic points by pre-setting deg=-1
    while ~orb(phi(a))
        a=phi(a);
        orb(a)=q;
        ord(a)=ord(psi(a))+1;
        deg(a)=-1;
    end
end
end % function orbits
%-----------------------------------------------------
function showmethemaps(SM)
  % showmethemaps(SM)
  % This is an internal auxiliar function to picture the endomaps from the matrix SM which is interpreted as:
  % SM(i,:)=[ Q(i), i, endomap(i,:)] 
  % where Q(i) is the code color of the endomap i
 
  [nn,n]=size(SM); n=n-2;
g=exp(linspace(0,2*pi*(n-1)/n,n)*1i);
na=round(sqrt(nn)); nb=ceil(nn/na);
if nn>64, error('too big to plot; SM must have less than 64 rows'), end
for i=1:nn,
  subplot(na,nb,i),
  dg=g(SM(i,3:n+2))-g;
  plot(g,'.'), hold on,
  quiver(real(g),imag(g),real(dg),imag(dg));
  axis off, hold off
  title(num2str(SM(i,1:2)))
end
disp('yes yes')
end % showmethemaps
%} % end of commented block, starts 484-210=274 lines above this one