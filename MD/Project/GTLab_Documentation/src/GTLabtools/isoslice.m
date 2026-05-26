function [a0,a1,ilevel,z0,z1,t0,t1]=isoslice(dom,cod,color,xval)
% [a0,a1,i,z0,z1]=isoslice(dom,cod,color,xval)
% 16Oct2024 -- by martins.ferreira@ipleiria.pt
% given a graph in the complex plane:
% dom,cod:A-->Complex (complex coordinates of each vertex)
% color:A-->Naturals (identifying faces by having the same color, eventhough they may have different connected components)
% xval:I-->Reals (a real value as refence to the slicing)
% It returns a sequence of tuples (a0,a1,i,z0,z1) where each tuple is representing an arrow from z0 to z1 corresponding to the xval indexed by i (in I) and anchored in the two edges a0 and a1 (in A) such that:
% a) real(z0)=real(z1)=xval(i)
% b) there exist t0 and t1 such that
% z0=dom(a0)+t0*(cod(a0)-dom(a0))
% z1=dom(a1)+t1*(cod(a1)-dom(a1))
% the edge a0 is oriented from right to left whereas a1 is from left to right
% the vector (z0,z1) in the complex plane is completely inside the face of color (color(a0)==color(a1)) containing edges a0 and a1.
% Example of application at the end of file

% forcing column vectors
d=dom(:); c=cod(:); x=xval(:); q=color(:);

% defining  A-by-I matrices
dI=repmat(d,1,numel(x));
cI=repmat(c,1,numel(x));
Ax=repmat(x',numel(d),1);

% filtering A-by-I matrices
up=(real(dI)<=Ax & Ax<real(cI));
down=(real(cI)<=Ax & Ax<real(dI));
[aup,iup]=find(up);
[adown,idown]=find(down);
tup=(x(iup)-real(d(aup)))./real(c(aup)-d(aup));
tdown=(x(idown)-real(d(adown)))./real(c(adown)-d(adown));
zup=d(aup)+tup.*(c(aup)-d(aup));
zdown=d(adown)+tdown.*(c(adown)-d(adown));
% checking that abs(real(zup)-r0(iup))<eps
all(abs(real(zup)-x(iup))<eps)
all(abs(real(zdown)-x(idown))<eps)

% the sorting and chaining process
upsort=sortrows([q(aup) iup imag(zup) aup]);
downsort=sortrows([q(adown) idown imag(zdown) adown]);
% checking that faces and levels are consistent in upsort and downsort
all(upsort(:,1:2)==downsort(:,1:2))
domcod=[downsort, upsort(:,3:4)];

a0=domcod(:,4);
a1=domcod(:,6);
ilevel=domcod(:,2);
z0=x(ilevel)+1i*domcod(:,3);
z1=x(ilevel)+1i*domcod(:,5);

t0=(z0-dom(a0))./(cod(a0)-dom(a0));
t1=(z1-dom(a1))./(cod(a1)-dom(a1)); % tup and tdown could be reused instead of computing t0 and t1

% end of function

%{ 
%Example of application
Draw=1;
% Testing example the letter A
[g,phi]=linkfigs('A',1+2i)
% adding the middle inside A
g=[g ; -.5+.5i; 1i; -1+1i]
phi=[phi 11 9 10 ]'
% Adding a hat to A in order to have two different faces
g=[g ; [-.5+.5i; 1i; -1+1i]-1i]
phi=[phi; [10 11 9]'+3]
if Draw
h=linkdraw(g,phi,1:14,.5,0);
set(h,"showarrowhead",'off')
pause
end

%g=g*1i; % rotated 90 degrees
g=g.*exp(1i*45/360*2*pi); % rotate 45
dom=g; cod=g(phi);
xval=linspace(-2,1,15);
faces=[1 1 2];
color=faces(orbits(phi))

% calling isoslice
[a0,a1,ilevel,z0,z1]=isoslice(dom,cod,color,xval);

% checking with linkdraw
if Draw
h=linkdraw(g,phi,[],0,0);
set(h,"showarrowhead",'off')
hold on
dir=z1-z0;
quiver(real(z0),imag(z0),real(dir),imag(dir),0)
text(real(z0+.5*dir),imag(z0+.5*dir),num2str((1:numel(ilevel))'))
hold off
pause
end

% reordering taking into account the connected components of phi
conn=orbits(phi)
[~, sortedbyconncomp]=sortrows([q(a0) conn(a0) ilevel z0 z1]);
% re-re-checking with linkdraw
if Draw
h=linkdraw(g,phi,[],0,0);
set(h,"showarrowhead",'off')
hold on
dir=z1-z0; dir_=dir(sortedbyconncomp);
z0_=z0(sortedbyconncomp); 
quiver(real(z0_),imag(z0_),real(dir_),imag(dir_),0)
text(real(z0_+.3*dir_),imag(z0_+.3*dir_),num2str((1:numel(ilevel))'))
hold off
pause
end

%}  
% end of file









