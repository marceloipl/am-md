% To evolve into a function constructing platonic solids via stereographic projection

Trial=3;  %
switch Trial
case 1 % Trial 1 (Thu, 21Nov2024; 14:15) 
% the thetrahedron
N=6; % number of edges on the base of the tetrahedron (seen as a pyramid)
% find the right value for the radius of the base circle
theta=@(z) acos(z) - 2 / N * pi * sqrt (1 - z .^ 2);
z=-fzero(theta,[-0.7,.9]);
r=sqrt(1-z^2);
pts=[r*exp((1:N)/N*2*pi*1i) , 0].';
phi_base=[N 1:N-1];
idA=(1:4*N)';
phi=zeros(N+3*N,1);
s=zeros(N+3*N,1);
phi(1:N)=[2:N 1];
s(1:N)=(1:N)';
for i=1:N
    phi(i+[N 2*N 3*N])=i+[ 2*N 3*N N];
    s(i+[N 2*N 3*N])=[ i N+1  phi_base(i)];
end
pos=pts(s);
V=plane2sphere(pos);
V(:,3)=-V(:,3);
linkdraw(V,phi,idA,.25,0)
case 2 % Trial 2 (Thu, 21Nov2024; 18:30)
N=6;
idA=(1:N)';
phi=[2:N 1]'; pos=exp(idA/N*2*pi*1i);

% from given N or given (pos,phi,idA) construct:
idA_=[idA; idA+N; idA+2*N; idA+3*N ];
phi_=[phi; idA+2*N; idA+3*N; idA+N ];
s_=[idA; phi; idA; repmat(N+1,N,1) ];
pos_=[pos-1/N*sum(pos); 0];
pos_=pos_/mean(abs(pos_(1:end-1)));

% visualize in 2d
subplot(1,2,1)
linkdraw(pos_(s_),phi_,idA_,.25,0);

% visualize with an offset
gh=linkoffset(idA_',phi_',pos_(s_).',.1).';
%linkdraw(gh,phi_,idA_,.25,0);


% visualize in 3d
subplot(1,2,2)
linkdraw(complex2sphere(pos_(s_)),phi_,idA_,.25,0);

case 3 % The cube [22Nov2024]
% given N do:
N=5;
phi=[2:N, 1]';
pos=exp((1:N)/N*2*pi*1i).';

% scale and recenter in the case when pos and phi are given arbitrarily, in which case we have to set
% N=numel(phi);
idA=(1:N)';
center=pos-sum(pos)/N;
pos=center/mean(abs(center));

% defining the cube out of a (generalized) base square
pos_=[1/2*pos; 2*pos];
phinv=zeros(numel(phi),1); phinv(phi)=idA;
idA_=[idA; idA+N; idA+2*N; idA+3*N; idA+4*N; idA+5*N];
phi_=[phi; phinv+N; idA+3*N; idA+4*N; idA+5*N; idA+2*N];
s_=[idA; phi+N; phi+N; phi; idA; idA+N];

% visualize in 2d and 3d
subplot(1,2,1)
h1=linkdraw(pos_(s_),phi_,s_,0,0);
set(h1,'showarrowhead','off')
% visualize in 3d
subplot(1,2,2)
h2=linkdraw(complex2sphere(pos_(s_)),phi_,s_,0,0);
set(h2,'showarrowhead','off')
end % Trials