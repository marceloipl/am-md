disp("Conjugation classes from p. 47 of Mat4CompSci.tex");
%% endomapsclass234students.m (last revised 5Dez2023)
%% Script to compute the conjugation classes of endomaps on a set with less than 5 elements
 
% Note that this code is only good for n=2,3,4
 
clf, clear, clc
n=3, nn=n^n, U=1:n; T=1:nn;
 
M=fillzeros(zeros(1,n),1:n)';
P=perms(1:n);
[x,y]=ndgrid(1:prod(1:n),1:n);
Pbar(sub2ind(size(P),x,P))=y;
Pbar=reshape(Pbar, size(P)); % the inverses of each P(i,:)
 
[~,S]=ismember(P,M,'rows');
[~,Sbar]=ismember(Pbar,M,'rows');
 
% computing ST as ST(s,t)=M(s,M(t,u)) linearized with s in S, t in T, and u in U
[s t u]=ndgrid(S,T,U); np=log(sqrt(exp(2*gamma(n+2))))/(n+1);
[~,ST]=ismember(reshape(M(sub2ind([nn,n],s,M(sub2ind([nn,n],t,u)))),[np*nn,n]),M,'rows');
ST=reshape(ST,[np,nn]);
 
% computing TSbar as TSbar(t,sbar)=M(t,M(sbar,u)) linearized with sbar in Sbar, t in T, and u in U
[tbar sbar ubar]=ndgrid(T,Sbar,U);
[~,TSbar]=ismember(reshape(M(sub2ind([nn,n],tbar,M(sub2ind([nn,n],sbar,ubar)))),[np*nn,n]),M,'rows');
TSbar=reshape(TSbar,[nn,np]);
 
Q=zeros(n,1);
flagcontinue=true; i=1;
while flagcontinue
[x,t]=ndgrid(1:np,i);
ST_i=ST(sub2ind(size(ST),x,t));
STSbar_i=TSbar(sub2ind(size(TSbar),ST_i,x));
        Q(STSbar_i)=i;
    fi=find(~Q,1); %  finds the first non-zero element in Q
    if isempty(fi)
        flagcontinue=false;
      elseif  fi==i
        error(['Entering in a cyclic loop at i = ' num2str(i)])
      else
        i=fi;
    end
end
 
 
% display the solutions
SM=sortrows([Q,(1:nn)',M]);
disp('[ Q(i), i, endomap(i,:)] ')
disp(SM)
 
% Picturing the results
if n<=3
  figure(1,"Name",['All cases n=' num2str(n)]), clf
showmethemaps(SM)
  figure(2,"Name",['Classes n=' num2str(n)]), clf
  [u,r,q]=unique(SM(:,1));
showmethemaps(SM(r,:))
elseif n==4
b1to4=[1 64 128 192 256];  %
b1=(1:64); b2=(65:128); b3=(129:192); b4=(193:256);
figure(1,"Name",['Case n=' num2str(n) '; (1:64)']), showmethemaps(SM(b1,:));
figure(2,"Name",['Case n=' num2str(n) '; (65:128)']), showmethemaps(SM(b2,:));
figure(3,"Name",['Case n=' num2str(n) '; (129:192)']), showmethemaps(SM(b3,:));
figure(4,"Name",['Case n=' num2str(n) '; (193:256)']), showmethemaps(SM(b4,:));
  figure(5,"Name",['Classes n=' num2str(n)]), clf
  [u,r,q]=unique(SM(:,1));
showmethemaps(SM(r,:))
else
disp('Too big to show the maps, use the function showmethemaps with smaller blocks instead')
end
 
% Checking that  isequal(ST,Tab(S,:)) and isequal(TSbar,Tab(:,Sbar)) where Tab is the full table
if n<4
Tab=zeros(n);
for i=1:nn
  for j=1:nn
    [~,pos]=ismember(M(i,M(j,:)),M,'rows');
    Tab(i,j)=pos;
  end
end
%Tab
isequal(ST,Tab(S,:))
isequal(TSbar,Tab(:,Sbar))
end
 
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