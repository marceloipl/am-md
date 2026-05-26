function [R,x,y]=hassediv(n,varargin)
% [R,x,y]=hassediv(n)
% 8Apr2025
% Given n, returns the partial order iRj iff i divides j, where i,j are two divisors of n
% in addition it returns coordinates [x,y] for the corresponding Hasse diagram
%
% See also the paper:
% https://math.hawaii.edu/~ralph/Preprints/latdrawing.pdf
% [2] Freese, Ralph (2004), Automated lattice drawing, Concept Lattices, LNCS 2961, Springer-Verlag, 589–590.
% [1] Di Battista, G.; Tamassia, R. (1988), Algorithms for plane representation of acyclic digraphs, Theoretical Computer Science 61 (2–3): 175–178.
% references based on wikipedia
% As well as Matlab Exchange Files:
% Andreas Plattner (2025). hasse diagram (https://www.mathworks.com/matlabcentral/fileexchange/27803-hasse-diagram), MATLAB Central File Exchange. Retrieved April 29, 2025.
% 
% Examples:
% R=hassediv(1000)
% [R,x,y]=hassediv(2025)
% hassediv(0,R) % R is the adjacency matrix of a partial order

if nargin==1
D=find(mod(n,1:n)==0); nD=numel(D);
[i,j]=ndgrid(1:numel(D));
R=(mod(D(j),D(i))==0); % D(i) divides D(j) i.e. R(i,j)==1
else % use R as varargin{1}
    R=varargin{1}
    D=1:size(R,1); nD=numel(D);
end
s=sum(R)'-sum(R,2);
[~,~,q]=unique(s);
[qs,sq]=sort(q);
disp([s(sq) qs   mod((1:nD)'-1,qs)+1 (1:nD)'])
y(sq)=s(sq);
v=hist(qs,1:max(qs))';
x(sq)=(mod((1:nD)'-1,v(qs))+1) - v(qs)/2  + 0.5*sq/nD;

y=(y-min(y))/(max(y)-min(y))
x=(x-min(x))/(max(x)-min(x))

% checking that R is reflexive, transitive and antisymmetric
isequal(ones(size(R,1),1),diag(R))
isequal(R,logical(R*R))
isequal(eye(size(R)),R.*R')

% in order to draw the lines in the hasse diagram
R2=(R^2==2);
[i,j]=find(R2);
du=x(j)-x(i);
dv=y(j)-y(i);
hq=quiver(x(i),y(i),du,dv,0);
 set(hq,'showarrowhead','off')
text(x,y,num2str(D(:))), axis([min(x) max(x) min(y) max(y)]), axis off


% in order to check if Rn is isomorphic to Rm
div=@(n) find(mod(n,1:n)==0);
Dn=div(20)
Dm=div(18)
nD=numel(Dn); mD=numel(Dm)
[in,jn]=ndgrid(1:nD);
[im,jm]=ndgrid(1:mD);
Rn=(mod(Dn(jn),Dn(in))==0);
Rm=(mod(Dm(jm),Dm(im))==0);

% f:Dn-->Dm, 
f=[1 3 5 2 4 6 ]
% g:Bm-->Dn, g=f^{-1}
g(f)=1:6;

% check if f is monotone
[i,j]=find(Rn);
all(Rm(sub2ind(size(Rm),f(i),f(j))))

% check if g is monotone
[i,j]=find(Rm);
all(Rn(sub2ind(size(Rn),g(i),g(j))))







% ------------------------------------------
% several initial experiments
%{
D=find(mod(n,1:n)==0); nD=numel(D);
[i,j]=ndgrid(1:numel(D));
R=(mod(D(j),D(i))==0); % D(i) divides D(j) i.e. R(i,j)==1
s=sum(R)'-sum(R,2);
%[u,r,q]=unique(s);
%v=hist(q,1:max(q))';
%{
tosort=[(s+nD), min(sum(R)',round(nD/2)), max(round(nD/2),nD-sum(R,2)), D(:), (1:nD)']
sorted=sortrows(tosort)
a=.75; b=.45; c=.35;
rms=round(1/3*sorted(:,1:3)*[.5*a 2*b c]')
rms=rms(sorted(:,5)); % make rms compatible with D
[~,~,q]=unique(rms);
v=hist(q,1:max(q))';
x_=mod((1:nD)',v(q))
x=x_-(v(q)-1)/2;
y=q-1;
%text(x,y,num2str(D(:)))
%text(x,y,num2str(sorted(:,4))), axis([min(x) max(x) 0 max(y)])
%}

% simplifying some aspects
[u,r,q]=unique(s);
[qs,sq]=sort(q);
disp([s(sq) qs   mod((1:nD)'-1,qs)+1 (1:nD)'])
y(sq)=s(sq);
%x(sq)=(mod((1:nD)'-1,qs)+1);
%x(sq)=log(mod((1:nD)'-1,qs)+1)
v=hist(qs,1:max(qs))';
%x(sq)=(mod((1:nD)'-1,qs)+1)-v(qs)/2; %-.5*nD./sq;
x(sq)=(mod((1:nD)'-1,v(qs))+1) - v(qs)/2  + 0.5*sq/nD;

y=(y-min(y))/(max(y)-min(y))
x=(x-min(x))/(max(x)-min(x))


% checking that R is reflexive, transitive and antisymmetric
isequal(ones(size(R,1),1),diag(R))
isequal(R,logical(R*R))
isequal(eye(size(R)),R.*R')

% in order to draw the lines in the hasse diagram
S=(R^2==2);
[i,j]=find(S);
du=x(j)-x(i);
dv=y(j)-y(i);
hq=quiver(x(i),y(i),du,dv,0);
 set(hq,'showarrowhead','off')
text(x,y,num2str(D(:))), axis([min(x) max(x) min(y) max(y)]), axis off
disp('end')
%}