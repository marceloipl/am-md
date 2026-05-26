function showmethemaps(SM)
  % showmethemaps(SM)
  % This is an internal auxiliar function to picture the endomaps from the matrix SM which is interpreted as:
  % SM(i,:)=[ Q(i), i, endomap(i,:)] 
  % where Q(i) is the code color of the endomap i
 
  [nn,n]=size(SM); n=n-2;
g=exp(linspace(0,2*pi*(n-1)/n,n)*1i);
na=round(sqrt(nn)); nb=ceil(nn/na);
if nn>64, error('too big to plot; SM must have less than 64 rows'), end
for i=1:nn,
  subplot(na,nb,i),
  dg=g(SM(i,3:n+2))-g;
  plot(g,'.'), hold on,
  quiver(real(g),imag(g),real(dg),imag(dg));
  axis off, hold off
  title(num2str(SM(i,1:2)))
end
disp('yes yes')
end