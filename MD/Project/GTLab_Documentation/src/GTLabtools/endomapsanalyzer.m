function [Rstar,p,q]=endomapsanalyzer(k)
%  [Rstar,p,q]=endomapsanalyzer(k)
% 30Abr2025; 10:00 --> see endomapsanalyzer30Abr
% 
% Possible solution to Class Activity proposed in Guide4Classes.pdf version 3

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
[a i]=ndgrid(1:nA,1:numel(d))
b=s2i(composerows(i2s(d(i)),
composerows(i2s(a),i2s(c(i)))));
S=logical(sparse(a,b,1,nA,nA))
%[a,b]=find(S);

[Rstar,p,q]=analyzerelation(S);
% analyzing the smallest equivalence relation containing S
linkdraw(homology(q),q)
pause
SM=sortrows([q, (1:nA)' i2s(1:nA)]);
showmethemaps(SM)
pause
% analyzing the partial order induced by S
u=unique(p)
Rstar(u,u)
analyzerelation(Rstar(u,u));
[x,y]=hasse(Rstar(u,u),u);


end % endomaps main function
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

fprintf('is reflexive? %d \n',all(diag(R)))
fprintf('is transitive? %d \n',all(all(logical(R*R)<=R)))
fprintf('is symmetric? %d \n',isequal(R,R'))
fprintf('is anti-symmetric? %d \n',all(all((R&R')<=eye(size(R)))))

Rtran=transitiveclosure(R);
Rstar=Rtran | eye(size(R)) % reflexive and transitive closure
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

y=(y-min(y))/(max(y)-min(y))
x=(x-min(x))/(max(x)-min(x))

% in order to draw the lines in the hasse diagram
R2=(R^2==2);
[i,j]=find(R2);
du=x(j)-x(i);
dv=y(j)-y(i);
disp([x(i),y(i),du,dv])
hq=quiver(x(i),y(i),du,dv,0);
 set(hq,'showarrowhead','off')
text(x,y,num2str(D(:))), axis([min(x) max(x) min(y) max(y)]), axis off
end % hasse