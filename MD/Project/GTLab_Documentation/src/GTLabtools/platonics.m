function [g,varphi,s,vstereo,vplain]=platonics(plato,pos,phi)
% [g,varphi,s]=platonics(type,pos,phi)
%
% This function is fully doccumented on Scripta-Ingenia, Number 13, December, 2024.
% plato is a number 1:5 where 1 is tetrahedron, 2 is cube, 3 is octahedron, 4 is dodecahedron and 5 is icosahedron (not yet implemented)
% (pos,phi) is a complex link encoding a base curve
% if phi is not given then cyclic order is assumed
% if pos is a number then a regular polygon with N=round(abs(pos)) sides is created as the base curve (pos,phi)
% if pos is not given then N=3, 4, 4, 5,3, respectively for plato=1,2,3,4,5

% controling user inputs
switch nargin()
case 1
 N=0; pos=[]; phi=[];
case 2
phi=[];
    if numel(pos)==1  % use pos(1) to set N
        N=round(abs(pos)); 
        pos=[];
    else
        N=numel(pos);
    end
case 3
N=numel(phi);
end

if any(strcmp(plato,{'tetra','tetrahedron','pyramid'})), plato=1; end
if any(strcmp(plato,{'cube','hexahedron','prism'})), plato=2; end
if any(strcmp(plato,{'octa','octahedron','double pyramid'})), plato=3; end
if any(strcmp(plato,{'dode','dodecahedron','pentagon'})), plato=4; end
if any(strcmp(plato,{'ico','icosahedron','star of triangles'})), plato=5; end

% Adjust N if necessary
Nref=[3 4 4 5 3];
if N==0 
    N=Nref(plato);
end
% N is updated, construct phi,idA, pos if necessary
if numel(pos)==0 % need to construct phi, idA, pos
phi=[2:N 1]'; idA=(1:N)'; pos=exp(idA/N*2*pi*1i); 
elseif numel(phi)==0 % need to construct phi, idA
    phi=[2:N 1]'; idA=(1:N)'; 
else % need to construct idA
    idA=(1:N)';
end
% pos, idA, phi, N, % are now given


switch plato

case 1 %{1, 'tetra'}
% from given N, pos, phi, idA  construct:
idA_=[idA; idA+N; idA+2*N; idA+3*N ];
phi_=[phi; idA+2*N; idA+3*N; idA+N ];
s_=[idA; phi; idA; repmat(N+1,N,1) ];
pos_=[pos-1/N*sum(pos); 0];
pos_=pos_/mean(abs(pos_(1:end-1)));
% output
g=pos_; varphi=phi_; s=s_;
vstereo=complex2sphere(g); % stereographic projection
vplain=[real(g) imag(g) [zeros(N,1); -1]];

case 2 %{2, 'cube', 'hexa'}
% from given N, pos, phi, idA  do:
% defining the cube out of a (generalized) base square
phinv=zeros(numel(phi),1); phinv(phi)=idA;
idA_=[idA; idA+N; idA+2*N; idA+3*N; idA+4*N; idA+5*N];
phi_=[phi; phinv+N; idA+3*N; idA+4*N; idA+5*N; idA+2*N];
s_=[idA; phi+N; phi+N; phi; idA; idA+N];
% scale and recenter in the case when pos and phi are given arbitrarily, in which case we have to set
centered=pos-sum(pos)/N;
scaled=centered/mean(abs(centered));
upper=[1/2*scaled];
lower=[2*scaled];
pos_=[upper; lower];
% output
g=pos_; varphi=phi_; s=s_;
vstereo=complex2sphere(g); % stereographic projection
vplain=[real(g) imag(g) [ones(N,1); -1*ones(N,1)]];

case 3 %{3, 'octa'}
% from given N, pos, phi, idA  construct:
idA_=[idA; idA+N; idA+2*N; idA+3*N; idA+4*N; idA+5*N ];
phi_=[idA+4*N; idA+2*N; idA+3*N; idA+N; idA+5*N; idA ];
s_=[idA; phi; idA; repmat(N+1,N,1); phi; repmat(N+2,N,1) ];
upper=N+2;
middle=idA;
lower=N+1;
centered=pos-sum(pos)/N; % the origin is the center of mass
scaled=centered/mean(abs(centered)); % the average radius is 1
pos_=[scaled; 0; NaN];  % 0=South Pole, NaN=North Pole
% output
g=pos_; varphi=phi_; s=s_;
vstereo=complex2sphere(g); % stereographic projection
vplain=[[real(scaled) imag(scaled) zeros(N,1)]; 0 0 -1; 0 0 1];

case 4 %{4, 'dodecahedron', 'penta'}
% from given N, pos, phi, idA  construct:
phinv=zeros(numel(phi),1); phinv(phi)=idA;
idA_=[idA; idA+N; idA+2*N; idA+3*N; idA+4*N; idA+5*N; ...
idA+6*N; idA+7*N; idA+8*N; idA+9*N; idA+10*N; idA+11*N];
phi_=[phi; idA+2*N; idA+3*N; idA+4*N; idA+5*N; idA+N; ...
idA+7*N; idA+8*N; idA+9*N; idA+10*N; idA+6*N; phinv+11*N];
s_=[idA; phi; idA; idA+N; idA+2*N; phi+N; idA+3*N; phi+3*N; phi+2*N; phi+N; idA+2*N; phi+3*N ];
upper=idA;
middleup=[idA+N; phi+N];
middledown=[idA+2*N; phi+2*N];
lower=[idA+3*N; phi+3*N];
centered=pos-sum(pos)/N; % the origin is the center of mass
scaled=centered/mean(abs(centered)); % the average radius is 1
pos_=[1/2*scaled; 2/3*scaled; 3/2*scaled; 2*scaled]; 
% output
g=pos_; varphi=phi_; s=s_;
vstereo=complex2sphere(g); % stereographic projection
plain=[scaled; scaled; scaled*exp(1i*2*pi/N); scaled*exp(1i*2*pi/N)]; 
vplain=[real(plain) imag(plain) [repmat(1,N,1); repmat(2/3,N,1); repmat(1/3,N,1); repmat(0,N,1)]];

case 5 %{5, 'ico', 'icosahedron'}
disp('Not implemented')
g=pos; varphi=phi; s=idA; vstereo=[]; vplain=[];

end


%{
%List of Examples
% Tetrahedron
[g,phi,s,v1,v2]=platonics(1,3);
subplot(4,3,1), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,2), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,3), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(1,4);
subplot(4,3,4), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,5), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,6), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(1,5);
subplot(4,3,7), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,8), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,9), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(1,6);
subplot(4,3,10), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,11), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,12), linkdraw(v2(s,:),phi,0), axis off

% The cube
[g,phi,s,v1,v2]=platonics(2,3);
subplot(4,3,1), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,2), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,3), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(2,4);
subplot(4,3,4), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,5), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,6), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(2,5);
subplot(4,3,7), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,8), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,9), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(2,6);
subplot(4,3,10), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,11), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,12), linkdraw(v2(s,:),phi,0), axis off

% The octahedron
[g,phi,s,v1,v2]=platonics(3,3);
subplot(4,3,1), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,2), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,3), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(3,4);
subplot(4,3,4), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,5), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,6), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(3,5);
subplot(4,3,7), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,8), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,9), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(3,6);
subplot(4,3,10), linkdraw(g(s,:),phi,3,.25), axis off
subplot(4,3,11), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,12), linkdraw(v2(s,:),phi,0), axis off

% The Dodecahedron
[g,phi,s,v1,v2]=platonics(4,3);
subplot(4,3,1), linkdraw(g(s,:),phi,0,.25), axis off
subplot(4,3,2), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,3), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(4,4);
subplot(4,3,4), linkdraw(g(s,:),phi,0,.25), axis off
subplot(4,3,5), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,6), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(4,5);
subplot(4,3,7), linkdraw(g(s,:),phi,0,.25), axis off
subplot(4,3,8), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,9), linkdraw(v2(s,:),phi,0), axis off

[g,phi,s,v1,v2]=platonics(4,6);
subplot(4,3,10), linkdraw(g(s,:),phi,0,.25), axis off
subplot(4,3,11), linkdraw(v1(s,:),phi,0), axis off
subplot(4,3,12), linkdraw(v2(s,:),phi,0), axis off

%}