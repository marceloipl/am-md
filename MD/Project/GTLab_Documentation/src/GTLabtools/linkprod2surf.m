function [Q,V,turn]=linkprod2surf(d,alpha,G,Phi,F)
% [Q,V]=linkprod2surf(phi1,g1mphi2,g2)
%
% Given a family of n complex link structures (indexed by u=(1:n) all of which
%  having m points, indexed by v=(1:m)), represented as (G(u,:),Phi(u,:)),
% together with a scalar function F(u,:) mapping each indexing point (u,v) in
%  the canonical euclidean 3-space with coordinates
% (x,y,z):=(F(u,v),real(G(u,v)),imag(G(u,v)))
% this means that each (G(u,:),Phi(u,:)) is considered as a non-planar curve with
% F(u,:) measuring at each point the difference from  planar, in the sense that
% when F(u,v)=0 then the point (u,v) lies on the complex plane which is being
% identified with the plane yz.
% The input (d,alpha) is also a complex link but it is interpreted as a planar
% curve in the xy plane which gives the direction of slicing, whereas each
% (G(u,:),Phi(u,;)) is the contour slice at level u. The size of the link (d,alpha)
% is m and the formula to place each point (u,v) in 3-space corresponding to its
% position on the slice at level u is:
% V(u,v)=[real(P(u,v)), imag(P(u,v)) imag(G(u,v))];
% where P(u,v) is:
% P(u,v)=d(u)+F(u,v)*turn(u)+1i*real(G(u,v)).*turn(u);
% where turn is the turning angle on the xy plane at each point u, that is,
% from d(u) to d(alpha(u)) to d(alpha(alpha(u)))
% It outputs a matrix V whose lines are the coordinates of each vertex, linearly
% indexed by (u,v) and  a matrix with quadrangles Q
%  in which each line represents a quadrangle
% connecting the four vertices in a cyclic order:
% Q(u,v)=[(u,v), (u,Phi(u,v)), (alpha(u),Phi(u,v)), (alpha(u),v)];
%{
% Simple example:
## n=7; m=2*n;
## tn=linspace(0,1,n);
## tm=linspace(0,1,m);
## d=exp(pi*tn*1i); alpha=[2:n,n];
## g=exp(2*pi*tm+1i); phi=[2,m,1];
## G=(1+.5*sin(n/4*pi*tn))*g(:).';
## Phi=repmat(phi,n,1);
## [Q,V]=linkprod2surf(d,alpha,G,Phi,zeros(n,m));
## tetramesh(Q,V)
%}

%[orb,ord]=orbits(alpha);
dir=d(alpha)-d; % direction vectors for the direction curve
dir=dir./(abs(dir));
dirac=dir(alpha)./dir; %dir and acceleration
turn=dir.*sqrt(dirac); % turning direction
turn=circshift(turn,1);
turn(isnan(turn))=0;
%disp(abs(turn))
%turn=circshift(turn,1);
%turn=cumprod([dir(1),sqrt(dirac(1:end-1))]);
%turn(turn==0)=turn(find(turn==0)-1); % turn(1) is assumed to be non-zero
%turn(turn==0)=interp1(find(turn~=0),turn(turn~=0),find(turn==0),'previous'); % replace zero values by interpolating neighbours
%turn=cumprod([dir(ord(1)),sqrt(dirac(ord(1:end-1)))]);
n=numel(alpha); m=numel(Phi)/n;
[u,v]=ndgrid(1:n,1:m);
q1=sub2ind([n,m],u,v);
q2=sub2ind([n,m],u,Phi);
q3=sub2ind([n,m],alpha(u),Phi);
q4=sub2ind([n,m],alpha(u),v);
Q=[q1(:),q2(:),q3(:),q4(:)];
P=d(u)+F.*turn(u)+1i*real(G).*turn(u);
x=real(P); y=imag(P); z=imag(G).*abs(turn(u));
V=[x(:),y(:),z(:)];
end

%{
% ------------------------Examples
disp('Example 1' )
 n=25; m=15;
 tn=linspace(0,1,n);
 tm=linspace(0,1,m);
 d=2*exp(1.5*pi*tn*1i); alpha=[2:n,n];
 g=exp(2*pi*tm*1i); phi=[2:m,1];
 s=1+0.2*sin(10*pi*tn);
 G=s'*g;
 Phi=repmat(phi,n,1);
 [Q,V]=linkprod2surf(d,alpha,G,Phi,zeros(n,m));
 tetramesh(Q,V), xlabel('x'), ylabel('y')
disp('Example 2')
d=[1,2,3,4,5,6]; alpha=[2 3 4 5 6 6];
phi2=[2:8,1]; phi3=phi2; phi3([2,6])=[7,3];
Phi=[phi2; phi2; phi3; phi3; phi2; phi2];
g2=[0 1 2 3 3+1i 2+1i 1+1i 1i];
g3=[0 1+.35i 2+.35i 3 3+1i 2+.65i 1+.65i 1i];
G=[g2;g2;g3;g3;g2;g2];
 [Q,V]=linkprod2surf(d,alpha,G,Phi,zeros(size(Phi)));
 tetramesh(Q,V), xlabel('x'), ylabel('y')
%}