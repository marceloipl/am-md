function [d,c,v]=hasse(R,varargin)
% [x,y]=hasse(R,D)
% 03May2025 
% 05Jun2025 a bug fixed: if the diagram is empty do not execute quiver
%
% given the adjacency matrix R of a partial order (a reflexive, transitive and antisymmetric relation) with elements D, return the x and y coordinates to draw the respective hasse diagram
% v=x+iy is the complex coordinate of each in D=varargin{1},
% [d,c] are the edges of the Hasse diagram,
% each edge i has coordinates v(d(i))-->v(c(i))
%
% Example:
% div=@(n) find(mod(n,1:n)==0)
% D=div(1000), [i,j]=ndgrid(1:numel(D)); R=(mod(D(j),D(i))==0);
% [x,y]=hasse(R,D)
% 
% interesting examples for hasse(n),
% n=1000,1008,5040,360,60

version='03May2025';
%version='30Abr2025'

if isequal(size(R),[1 1])
    div=@(n) find(mod(n,1:n)==0);
 n=R; % if R is of size 1-by-1 interpret it as R(i,j) iff i divides j
 D=div(n), [i,j]=ndgrid(1:numel(D)); R=(mod(D(j),D(i))==0);
elseif nargin==1, 
    D=1:size(R,1);
else
    D=varargin{1};
end

switch version

case '03May2025'

nD=numel(D);
if ~isequal(size(R),[nD nD]), error('R must be nD-by-nD, with nD=numel(D)'), end
if ~all(all(logical(R*R)<=R)) | ~all(diag(R)) | ~all(all((R&R')<=eye(size(R))))
    error('R is not a partial order, i.e. a reflexive, transitive and antisymmetric relation')
end

y=ones(nD,1);
R1=R*R==2; % paths of length 1, counting with loops
[d,c]=find(R1);
y(c)=y(c)+1;
Rn=R1*R1; % Get the number of paths of length n
while nnz(Rn)>0
    [i,j]=find(Rn);
    y(j)=y(j)+1;
    Rn=Rn*R1;
end
w=sparse(y,1,1);
[ys,sy]=sort(y); %ys=ysorted, sy=sorted indexes of y
cw=cumsum(w);
x(sy)=(1:nD)'-cw(ys)+w(ys)/2;
x=x(:);
v=x+1i*y;
%{
w=hist(y,1:max(y))';  % width of each layer
out=acumarray(y)
%S= sortrows([w(y) y sum(R,1)' -sum(R,2) D(:)  (1:nD)']);
S= sortrows([w(y) y  D(:)  (1:nD)'])
w_=S(:,1); id_=S(:,end);
x(id_)=mod((1:nD)'-1,w_)+1-(w_-1)/2
%x(S(:,3))=mod(S(:,3)-1,w(S(:,1)))+1 - w(S(:,1))/2  + 0.5*S(:,3)/nD;
%x(S(:,3))=mod(S(:,3)-1,w(S(:,1)))+1 - w(S(:,1))/2;
x=x(:);

% alternatively
%x=mod((1:nD)'-1,w(y))+1-(w(y)-1)/2
%}

% in order to draw the lines in the hasse diagram
du=x(c)-x(d);
dv=y(c)-y(d);
% disp('input to quiver'), disp([x(d),y(d),du,dv])
if ~isempty(d)
hq=quiver(x(d),y(d),du,dv,0);
 set(hq,'showarrowhead','off')
end
text(x,y,num2str(D(:))), axis([min(x) max(x) min(y)-eps max(y)+eps]), axis off


%------------------------------------------------------------------
case '30Abr2025'

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

end % switch
end % hasse