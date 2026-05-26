function theta=dlink(g,phi,s,nA)
% theta=dlink(g,phi,s,nA)
%
% 30Set2024 --- in progress (not yet complete) -- discontinued
% given the incomplete double link structure (g,phi,s,nA) where g:A-->Complex, phi:A-->A and s:A-->B where A=(1:nA) and B=(1:max(s)), where each a in A is interpreted as an arrow from s(a) to s(phi(a)) and g(a)=r(a)*exp(1i*t(a)) with r(a) the acceleration from a to phi(a) while t(a) in ]0,pi[ is the turning angle from a to phi(a), also refered as ampli-twist (r is the ampli while t is the twist)
% returns theta:A-->A such that theta(phi(theta(phi)))=1_A and s(theta)=s

% Examples
[V,T] = platonic_solids(2);  % dodecahedron
[g,phi,s]=tri2link(T,V)


try
idA=(1:nA)';
catch
idA=(1:numel(g))';
end

[~,theta_phi]=ismember([s(phi) s],[s,s(phi)],'rows');
theta(phi)=theta_phi; theta=theta(:);
[s_,sord]=orbits(theta); % the orbits of theta ate the vertices
[q,qord]=orbits(phi); % the orbits of phi are the faces
idV=(1:max(s_))';  % identity on vertices
idF=(1:max(q))';
%idA=(1:numel(g))';

Vord=sortrows([s_ sord, idA]);
Ford=sortrows([q qord, idA]);

g_=cumsum([0; exp(cumsum([0; log(g(Ford(1:end-1,3)))]))])

r=abs(g(Ford(:,3))); % ampli
t=angle(g(Vord(:,3))); % twist

defang=full(cumsum(sparse(1,s(phi(Vord(:,3))),pi-t)))';
refang=2*pi*Vord(:,1).*cumsum(pi-t)./defang(Vord(:,1));

A=full([sparse(q,idA,1); sparse(s(phi),idA,1); sparse(orbits(theta(phi)),idA,1/2)]);
% matrix A is such that if given f0 on vertices and faces, say fB and fC with B the index set of vertices and C the index set of faces, we find f distributed by edges so that f0 is recoverd from f
% example 
fB=[1+1i; zeros(6,1); -1-1i]; fC=0*idF;
f=A\[fB; fC; zeros(numel(idA)/2,1)];




% Mass matrix -- to be computed with similar tehcnique to refang

% The connection:
Civitta=exp(1i*(refang - refang(theta(phi))) + pi);
% We could now redefine W as
W=sparse(s,s(phi),abs(g).*cot(angle(g)/2).*Civitta);
W=0.5*(W+W'); % force symmetry and get the equivalent 1/2(cot(a)+cot(b)) -- however, it is not clear which effect will it have on the imaginary part
L=W-diag(sum(W));  %disp(full(L))



%{Examples
[v,f]=platonic(2); % a cube
[g,phi,s]=tri2link(f,v)

[V,T] = platonic_solids(5);  % dodecahedron
[g,phi,s]=tri2link(T,V)

%}