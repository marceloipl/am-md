function [nA,nV,nF,theta,M,L,divergence,rotational,grad]=dlinktopology(g,phi,s)
% [nA,nV,nF]=dlinktopology(g,phi,s)
%
% created: 21jun2024
% revised 20Set2024
%
% given [g,phi,s] as from tri2link, returns [nA,nV,nF,theta] with nA, nV, nF respectively the number of edges, vertices and faces whereas theta is the completion of (phi,s) into a double link (theta,phi,s,orb) with orb the orbits of phi while s can be interpreted as the orbits of theta
% the genus is: nV-nA/2+nF
% Examples at the end of file


nA=numel(phi); % number of edges

% preliminaries
id=@(x) (1:numel(x))';
phinv(phi)=id(phi);
w=-1./g(phinv); % g is such that if edge a in A has geometry a=r*exp(1i*t) then phi(a) has geometry a*g(a)
integralarea=@(a,b,c,d) 1/6*(b-a).*(d.^2+d.*c+c.^2);
% backward area, uses w(a)=1/g(phinv(a)) % 20jun but fixed here w(a)=-1/g(phinv(a))
backward=integralarea(0,angle(w),1,abs(w));
forward=integralarea(angle(g),pi,abs(g),1);
% backward area is used for the mass matrix M, whereas forward area is used in diverg

% The mass matrix, diverg and grad
M=sparse(s,1,backward);
divergence=@(X) sparse(s,1,forward.*real(X));  % X is a complex field indexed on edges 
rotational=@(X) sparse(s,1,forward.*imag(X));  % X is a complex field indexed on edges % to be tested
grad=@(u) 1i*( u(s).*g + u(s(phi)).*(-1-g) + u(s(phi(phi))) );   % u is a scalar filed indexed on vertices 

% completing the double-link structure: (theta,phi)
[~,theta]=ismember([s s(phinv)],[s s(phi)],'rows'); % see heat2dist16jun.m line 51
[orb,ord]=orbits(theta);
[~,alphaord]=sort(ord);  % alphaord is the permutation that sorts the vector ord, thus producing cyclic sequences around each vertex
% computing relative reference angles around each vertex
relrefang=360*round(3600*backward./M(s))/3600;   % working in degress with 1/3600 resolution --- this may be conveninet for cumsum to identify each beginning of a new vertex as a multiple of 2*pi
% computing the absolute reference angles around each vertex
refang=cumsum(relrefang(alphaord)*2*pi/360); % and returning to radians
% in principle the reference angles around vertex b in (1:numel(M)) are of the form x+k*2*pi with k=b, this is a consequence of cumsum but when we take exp(1i*refang) all will be well
% The connection:
Civitta=exp(1i*(refang - refang(theta(phi))) + pi);
% We could now redefine W as
W=sparse(s,s(phi),abs(g).*cot(angle(g)/2).*Civitta);
W=0.5*(W+W'); % force symmetry and get the equivalent 1/2(cot(a)+cot(b)) -- however, it is not clear which effect will it have on the imaginary part
L=W-diag(sum(W));  %disp(full(L))

q=coeq(id(phi),phi);
[u,r,q]=unique(q);
nV=numel(M);
nF=numel(r);  % gives the number of equivalence classes in the orbits of phi
%fprintf('The genus is %d', (nV-nA+nF-2)/2) --- is it?

% Some examples
%{
[T,V]=torusknot;
[g,phi,s]=tri2link(T,V);
%[L,M,grad,diverg]=laplacian(g,phi,s);
[nA,nV,nF,theta,M,L,diverg,rotational,grad]=dlinktopology(g,phi,s)
heat=@(u0,h) (diag(M)-h*L)\(diag(M)*u0(:));
dist=@(spikes,h)  sqrt(-4*h*log(heat(ismember(1:numel(M),spikes),h))); % Varadhan's formula
K=@(h) (12/h^2)*(1-M/pi);  % Gauss curvature
% in order to have the distance from Crane's heat method we observe the following steps
spike=1; h=1/sqrt(pi);
u=heat(ismember(1:numel(M),spike),h);
du=grad(u);
cranedist=L\diverg(du./abs(du));
disp('hello')
spikes=1:8; % heat source
d=dist(spikes,1);
trimesh(T,V(:,1),V(:,2),V(:,3),d), axis equal
hs=spikes; text(V(hs,1),V(hs,2),V(hs,3),num2str(round(d(hs)*10)/10))
%}