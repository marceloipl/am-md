function [V,Q]=links2surf(phi1,g1,phi2,g2,C0,C1)
% [V,Q]=links2surf(phi1,g1mphi2,g2)
%
% Given two complex link structures (phi1,g1) and (phi2,g2), in which the first 
% one is interpreted as a mother curve (the one that determines the center and position 
% of each copy of the generating curve) and the second one is interpreted as the geratrix. 
% The mother curve is identified with the plane xy while the genratrix is initially 
% identified with the yz plane which is rotated and placed at the places indicated by 
% the direction of the mother curve
% It outputs a matrix V whose lines are the coordinates of each vertex, indexed by 
% (i,j) with i ranging from 1 to numel(phi1) and j ranging from 1 to numel(phi2); 
% a matrix with quadrangles Q is also output in which each line represents a quadrangle 
% connecting the four vertices in a cyclic order
%
% C0 and C1 are complex vectors of the same size as g2 and serve as translation and rotation+spcaling to the curve g1.
% namely, for each j in the indexes of g2 we observe  

% Simple example:
% t=linspace(0,2*pi,7);
% g1=exp(t(1:end-1)*1i); phi1=[2:6,1];
% g2=[-.5 0 1i -0.5+1i]; phi2=[2:4,1];
% [V,Q]=links2surf(phi1,g1,phi2,g2);
% tetramesh(Q,V)

n=numel(phi1); m=numel(phi2);
dir=g1(phi1)-g1; % direction vectors of the mother curve
[i,j]=ndgrid(1:n,1:m);
a=sub2ind([n,m],i,j);
b=sub2ind([n,m],phi1(i),j);
c=sub2ind([n,m],phi1(i),phi2(j));
d=sub2ind([n,m],i,phi2(j));
Q=[a(:),b(:),c(:),d(:)];
th=angle(dir(i).*(dir(phi1(i))./dir(i)).^0.5); % angle of the vector tangent at g(phi1(i))
g2C01=C1(i).*g2(j)+C0(i);
x=real(g1(phi1(i)))+real(g2C01).*cos(th-pi/2);  % the angle th-pi/2 is ortogonal to the mother curve at each point g1(phi1)
y=imag(g1(phi1(i)))+real(g2C01).*sin(th-pi/2);
z=imag(g2C01);
V=[x(:),y(:),z(:)];
