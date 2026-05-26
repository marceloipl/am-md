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

% Examples
%{ 
z=1/4*[1,2,3+1i,3,4,4+6i,3+6i,3+2i,1+3i,2i,1i,1+1i,2+1i,2+2i,1+2i];
z=[z, sqrt(2)*exp(linspace(-pi,pi,16)*1i)];
% several cases can be analyzed:
% case 1: letter d!
zi=[0,1,1i]; wi=[0,1,1i];
w=mobious(z,zi,wi); plot(z), hold on, plot(w), hold off
% case 2: letter b?
zi=[0,1,1i]; wi=[0,-1,1i];
w=mobious(z,zi,wi); plot(z), hold on, plot(w), hold off
% case 3: letter q?
zi=[0,1,1i]; wi=[0,1,-1i];
w=mobious(z,zi,wi); plot(z), hold on, plot(w), hold off
% case 4: letter p?
zi=[0,1,1i]; wi=[0,-1,-1i];
w=mobious(z,zi,wi); plot(z), hold on, plot(w), hold off
% case 5: letter a?
zi=[3/4*[1i,1+1i], 1+3/2*1i]; wi=[3/4*[1i,1+1i], 1+3/4*1i];
w=mobious(z,zi,wi); plot(z), hold on, plot(w), hold off
% case 6: italic d?
zi=[0,1,1i]; wi=[0,1,1/4+1i];
w=mobious(z,zi,wi); plot(z), hold on, plot(w), hold off
%}

% More examples
%{
clf
[x,y]=ndgrid(-2:.5:2);
z=x+y*1i;
t=reshape(1:numel(x),size(x));
% case 1
subplot(1,2,1), cla, axis([-2 2 -2 2])
text(real(z),imag(z),num2str(t(:)),'fontsize',12),
w=mobious(z,[0 1 1i],[0 1 1i])
subplot(1,2,2), cla, axis([-2 2 -2 2])
text(real(w),imag(w),num2str(t(:)),'fontsize',12)
% case 2
subplot(1,2,1), cla, axis([-2 2 -2 2])
text(real(z),imag(z),num2str(t(:)),'fontsize',12),
w=mobious(z,[0 1 1i],[0 1.1 1i-.1])
subplot(1,2,2), cla, axis([-2 2 -2 2])
text(real(w),imag(w),num2str(t(:)),'fontsize',12)

% case 3
subplot(1,2,1), cla, 
%axis([-2 2 -2 2])
plot(z)
text(real(z),imag(z),num2str(t(:)),'fontsize',12),
w=mobious(z,[0 1 1i],[1 1+1i -1+1i])
subplot(1,2,2), cla,
wx=real(w(:)); wy=imag(w(:)); 
%axis([min(wx) max(wx) min(wy) max(wy)])
plot(w)
text(wx,wy,num2str(t(:)),'fontsize',12)
%}