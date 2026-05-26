function [X,Y]=sphere2plane(x,y,z)
% [X,Y]=plane2sphere(x,y,z)
% Using stereographic projection maps a point on the  unit sphere with coordinates (x,y,z) into the corresponding point on the plane with coordinates (X,Y) 

X=x./(1-z);
Y=y./(1-z);

%X(isnan(X))=Inf;
%Y(isnan(Y))=Inf;
