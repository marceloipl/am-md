function [V,x,y,z]=plane2sphere(X,Y)
% [x,y,z]=plane2sphere(X,Y)
% Using stereographic projection maps a point in the plane with coordinates (X,Y) into the corresponding point on the unit sphere with coordinates (x,y,z)

try Y; % if Y is not given assume X is complex
catch Y=imag(X); X=real(X); end

R=1+X.^2+Y.^2;
x=2*X./R;
y=2*Y./R;
z=(-1+X.^2+Y.^2)./R;

I=find(isnan(z));
x(I)=0;y(I)=0;z(I)=1;
V=[x y z];

%{ 
examples
[a b c]=plane2sphere(real(gA),imag(gA))

t=linspace(0,1,30);
V=[a b c]
plot3(a,b,c,'*')
hold on
for i=1:8
gp=slerp(V(i,:),V(phi(i),:),t)
plot3(gp(:,1),gp(:,2),gp(:,3))
end

% A study of the letter M and letter V
[phi,gM]=linkfigs('M',0);
[~,gV]=linkfigs('V',0);

t=linspace(0,1,6);
g=gM(:)*(1-t)+gV(:)*t;
for i=1:numel(t)
    linkdraw23(phi,g(:,i))
    pause
end

[a b c]=plane2sphere(real(gM),imag(gM));
V=[a(:) b(:) c(:)];
hold on
t=linspace(0,1,30);
for i=1:numel(phi)
   Vi=slerp(V(i,:),V(phi(i)),t);
   plot3(Vi(:,1),Vi(:,2),Vi(:,3))
   end
   hold off
%}
