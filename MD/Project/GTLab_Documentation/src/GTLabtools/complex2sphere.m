function v=complex2sphere(w)
% v=plane2sphere(w)
% see also plane2sphere
% Using stereographic projection maps a point in the complex plane with coordinates w=(X,Y) into the corresponding point on the unit sphere with coordinates v=(x,y,z)
% when w=NaN it returns v=[0 0 1]

X=real(w); Y=imag(w);
I=find(isnan(w));

R=1+X.^2+Y.^2;
x=2*X./R;
y=2*Y./R;
z=(-1+X.^2+Y.^2)./R;

x(I)=0;y(I)=0;z(I)=1;
v=[x y z];
