function [pos,phi,color]=linkpictures(Pict)
% [pos,phi,color]=linkpictures(Pict)
% returns the link strcuture for picture Pic

Draw=0;  % change to 0 to disable
Mobious=0; % change to 0 to disable
try
[pos,phi,color]=eval(Pict);
catch
switch Pict
case 'A'
[pos,phi]=linkfigA;
color=ones(numel(phi),1);
case '1'
[pos,phi]=linkfig1;
color=ones(numel(phi),1);
case '+'
[pos,phi]=linkfigplus;
color=ones(numel(phi),1);
case '++'
[pos,phi,color]=linkfigplusplus;
case 'Tree'
[pos,phi,color]=linkfigtree;
Ang=pi/12;
pos=pos*exp(1i*Ang);
case 'cyclicharmonics'
[pos,phi,color]=cyclicharmonics(4,1,100,3,1/3);
end
if Draw
    if Mobious
        wi=input('[0,1,1i]-->[w1,w2,w3]');
        try
        pos=mobious(pos,[0,1,1i],wi);
    end
    end
    h1=linkdraw(pos,phi,(1:numel(phi)),.5,0); axis equal
    yn=input('Draw trajectories? (1/0)');
    if yn
        xval=linspace(min(real(pos)),max(real(pos)),3*numel(phi));
        [a0,a1,ilevel,z0,z1]=isoslice(pos,pos(phi),color,xval);
        set(h1,"showarrowhead",'off')
hold on
dir=z1-z0;
h2=quiver(real(z0),imag(z0),real(dir),imag(dir),0)
set(h2,"showarrowhead",'off')
%text(real(z0+.5*dir),imag(z0+.5*dir),num2str((1:numel(ilevel))'))
hold off
end % if Draw
end % if yn
end % try
        
end % main function

%-------------------------------------------
function [pos,phi,color]=cyclicharmonics(p,q,N,a,e)
% [pos,phi]=cyclicharmonics(p,q,N,a,e)
%
% 21Mar2025
% implements cyclic-harmonics curves as in:
% https://mathcurve.com/courbes2d.gb/conchoidderosace/conchoidderosace.shtml
% and returns a complex-link structure [pos,phi]

n1=N; % using the same code as in torusknot.m

t1=linspace(1/n1,1,n1)*2*pi;

%if ~isequal(p,0), t1=t1/p; end

rho=a*(1+e*cos(p*t1));
x=rho.*cos(q*t1);
y=rho.*sin(q*t1); 

pos=x+1i*y;
phi=[2:N, 1];
color=ones(size(phi));
end % cyclicharmonoics

%-------------------------------------------
function [pos,phi]=linkfigA
a=pi/3;
b=1+1/cos(a);  %b=3
r1=[1 1 1 1 b 1/b b 1/b]; r2=[1 1 1]; 
t1=[a -a -a a pi-a a a pi-a ]; t2=[pi-a pi-a pi-a]; 
r=[r1 r2 ]';
t=[t1 t2 ]';
phi=[2:8 1 11 9 10]';
g=log(r)+1i*t;
dir(phi)=exp(cumsum(g));
pos(phi)=cumsum(dir);
pos(9:11)=pos(9:11)+(3+2*cos(a))/2+1i*2.5*sin(a);
pos=pos(:);
end % figA

%---------------------------------------------
function [pos,phi]=linkfig1
a=sqrt(3)/2;
pos=[0 2 1.5+a*1i 1.5+3.5*1i .5+3.5*1i (3.5-a)*1i .5+(3.5-a)*1i .5+a*1i].';
phi=[2:8,1]';
end % fig1

%----------------------------------------------
function [pos,phi]=linkfigplus
u=1; v=1i;
dir=[u v u v -u v -u -v -u -v u -v];
phi=[2:12,1]';
pos(phi)=cumsum(dir);
pos=pos.';
end % fig1

%----------------------------------------------
function [pos,phi,color]=linkfigplusplus
u=1; v=1i;
dir=[u v u v -u v -u -v -u -v u -v [u/2 v/2 -u/2 -v/2 ] ];
phi=[2:12,1 [14 15 16 13]]';
pos(phi)=cumsum(dir);
[conn,color]=linksignature([12 4],[1 2]);
pos(color==2)=pos(color==2)-1;
pos=pos.';
end % fig1

%-----------------------------------
function [pos,phi,color]=linkfigtree
u=1; v=1i;
Tree=[3*u v -u 4*v 3*u v u v -u v -u v -u v -u v -u v ];
pos=[cumsum(Tree )];
pos=[0 pos -conj(pos(end:-1:1))];
phi=[2:numel(pos),1]';
pos=pos(:);
color=ones(numel(pos),1);
end % tree


%-----------------------------------------------
function [conn,color]=linksignature(connsig,colorsig)
%connsig=[8 3 4];
%colorsig=[1 1 2];
n=[0 cumsum(connsig)];
conn=cumsum(ismember((1:n(end))',n(1:end-1)+1));
color=colorsig(conn);
end % auxilar funtion to return connected components and color faces from signatures

%------------------------------------------
function w = mobious(z,zi,wi)
% w = mobious(z,zi,wi) 
% computes w=f(z) where f(z) is the unique Mobious transformation f(z)=(az+b)/(cz+d) for which wi=f(zi) with zi=[z1, z2, z3] any three different complex numbers that are transformed into wi=[w1, w2, w3].

un=[1 1 1]';
zi=zi(:); wi=wi(:); % force column vectors
a=det([zi.*wi, wi, un]);
b=det([zi.*wi, zi, wi]);
c=det([zi, wi, un]);
d=det([zi.*wi, zi, un]);
w=(a*z+b)./(c*z+d);
end % mobious function

%---------------------------
function h=linkdraw(g,phi,labels,pos,arr)
% h=linkdraw24(phi,g,labels)
%
% New version 8Oct2024, the input g can be a complex vector as before but now it can also be a n-by-3 matrix with 3d coordinates
% plots the link structure with arrows using labels
% based on linkdraw24.m dated 29May2024 - the order of input agguments was changed: instead of phi,g now is g, phi, and if phi is not specified then it is assumed to be cyclic
% pos is the position of the labels, the default is pos=0.33
% arr is the arrow scaling, default is 0, or else not specified

%if iscomplex(g)  % in case g is complex 
if size(g,2) < 3 % is complex

try
dg=g(phi)-g;
catch
dg=g([2:end,1])-g;  % if phi is not given then assume cyclic order
end

try
h=quiver(real(g),imag(g),real(dg),imag(dg),arr);
%quiver(real(g),imag(g),real(dg),imag(dg),0); % to draw complete arrows
catch
h=quiver(real(g),imag(g),real(dg),imag(dg));
end

try
text(real(g+pos*dg),imag(g+pos*dg),num2str(labels(:)),"fontsize",12)
%text(real(g),imag(g),num2str(labels(:)))
catch
try
text(real(g+.33*dg),imag(g+.33*dg),num2str(labels(:)),"fontsize",12)
end
end

%---------The case g is a 3D vector
else
disp('using quiver3')

try
dg=g(phi,:)-g;
catch
dg=g([2:end,1],:)-g;  % if phi is not given then assume cyclic order
end
x=g(:,1); y=g(:,2); z=g(:,3);
dx=dg(:,1); dy=dg(:,2); dz=dg(:,3);
try
h=quiver3(x,y,z,dx,dy,dz,arr);
%  to draw complete arrows
catch
h=quiver3(x,y,z,dx,dy,dz); % if arr does not work simply ignore it
end

try
text(x+pos.*dx,y+pos.*dy,z+pos.*dz,num2str(labels(:)),"fontsize",12)
%text(real(g),imag(g),num2str(labels(:)))
catch
try
pos=.33;
text(x+pos.*dx,y+pos.*dy,z+pos.*dz,num2str(labels(:)),"fontsize",12)
end
end
end % end if else
end % function linkdraw