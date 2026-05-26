function w=sphere2complex(v)
% w=sphere2complex(v)
% see also [X,Y]=plane2sphere(x,y,z)
% Using stereographic projection maps a point on the  unit sphere with coordinates v=(x,y,z) into the corresponding point on the complex plane with coordinates w=X+1i*Y 
% if v=[0 0 1] it returns w=NaN

x=v(:,1); y=v(:,2); z=v(:,3);

X=x./(1-z);
Y=y./(1-z);

w=X+1i*Y;