function [T,V]=torusknot(p,q,n1,n2,radii)
%[T,V]=torusknot(p,q,n1,n2,radii)
% Tue, 4Mar2025 -- bug fixed by interchanging the roles of p and q; added parameters e and a that can be varied in the code
%
% older version dated: Sun, 9Jun2024
% [T,V]=torusknot(p,q,n1,n2,radii) 
% https://home.adelphi.edu/~stemkoski/knotgallery/
% https://mathcurve.com/courbes2d.gb/conchoidderosace/conchoidderosace.shtml
% given input numbers (p,q) returns a triangulation (T,V) for the (p,q)-torus knot such as for example the Trefoil (3,2) Torus Knot;
% the numbers n1 and n2 control the resolution n1-by-n2 of the quadrangular mesh, where n1 is the number of points in the main curve whereas n2 is the number of points in the cross-section, default values are n1=12 and n2=6
% radii is the radius of the cross-section which compares with 1 the radii of the main curve
%
% Simpler example:
% [T,V]=torusknot(0,1,3,3); trimesh(T,V(:,1),V(:,2),V(:,3)), axis equal
% More complex example:
% [T,V]=torusknot(3,2,30,13); trimesh(T,V(:,1),V(:,2),V(:,3)), axis equal
% [T,V]=torusknot(2,3,30,13); trimesh(T,V(:,1),V(:,2),V(:,3)), axis equal

a=3;
e=1/3; % e=1, e=3 
% see how to set different values at:
 % https://mathcurve.com/courbes2d.gb/conchoidderosace/conchoidderosace.shtml

switch nargin
case 0
p=0; q=1; n1=12; n2=6; radii=NaN;
case 1
q=1; n1=12; n2=6; radii=NaN;
case 2
n1=12; n2=6; radii=NaN;
case 3
n2=6; radii=NaN;
case 4
radii=NaN;
end

 t1=linspace(1/n1,1,n1)*2*pi;
 t2=linspace(1/n2,1,n2)*2*pi;

%if ~isequal(p,0), t1=t1/p; end

rho=a*(1+e*cos(p*t1));
 x=rho.*cos(q*t1);
 y=rho.*sin(q*t1); 
 z=q*sin(p*t1);
 
 g1=x+1i*y;
 phi1=[2:n1,1];
%r2=1-1./(log(n1).^n1);  % the ray of the cross-section
 
 if isnan(radii), radii=1./log(q+3); end % if radii is not specified then the formula 1./log(q+3) is used
 g2=radii*exp(1i*t2); 
 phi2=[2:n2,1]; 
 C0=1i*z; 
 C1=1+0*g1; 
 [V,Q]=links2surf(phi1,g1,phi2,g2,C0,C1); 
 T=[[Q(:,1) Q(:,2) Q(:,3)] ; [Q(:,1) Q(:,3) Q(:,4)]]; % convert quadrangulation into a triangulation
 
 end % function
%{ 
% older version dated 9Jun2024
 try p; q; n1; n2; catch p=1; q=0; n1=12; n2=6; end
 try n1; n2; catch n1=12; n2=6; end
 
 t1=linspace(1/n1,1,n1)*2*pi;
 t2=linspace(1/n2,1,n2)*2*pi;
 
 a=3; e=1/3; % if e=1 then cardioids are produced, see:
 % https://mathcurve.com/courbes2d.gb/conchoidderosace/conchoidderosace.shtml
 x=cos(p*t1).*(3+cos(q*t1));
 y=sin(p*t1).*(3+cos(q*t1)); 
 z=p*sin(q*t1);
 
 g1=x+1i*y;
 phi1=[2:n1,1];
 %r2=1-1./(log(n1).^n1);  % the ray of the cross-section
 r2=1./log(p+3);
 try r2=radii; end % if radii is not specified then r2 is used
 g2=r2*exp(1i*t2); 
 phi2=[2:n2,1]; 
 C0=1i*z; 
 C1=1+0*g1; 
 [V,Q]=links2surf(phi1,g1,phi2,g2,C0,C1); 
 T=[[Q(:,1) Q(:,2) Q(:,3)] ; [Q(:,1) Q(:,3) Q(:,4)]]; % convert quadrangulation into a triangulation
 
 end % function
%}
 
 %% Examples
%{
 [T,V]=torusknot(1,0,120,12);
 trimesh(T,V(:,1),V(:,2),V(:,3))
 text(V(:,1),V(:,2),V(:,3),num2str(round(D)))
 %}