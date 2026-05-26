function [phi,d_,c_,In,Out]=linkdigraph(d,c)
% [phi,d_,c_],In,Out,=linkdigraph(d,c)
% current version: 05May2025; created 29May2024
% support: martins.ferreira@ipleiria.pt
% 
% Given two vectors of the same size [d,c], interpreted as a digraph d,c:A-->B, returns a map phi:A-->A such that c=d(phi)
% In addition it returns In(b) and Out(b), the In and Out degree of each vertex b in B, as well as the new digraph [d_,c_] where d_ and c_ may contain extra edges so that the linking c_=d_(phi) is possible


version='05May2025';

switch version
case '05May2025'

%{
% Examples:
[d,c,v]=hasse(1008)
[phi,d_]=linkdigraph(d,c)
linkdraw(v(d_),phi,0)
linkdraw(v(d_),phi,3,0,d_)
%}

%d=[1 2 2 4 2 5], c=[2 3 4 3 2 4]
%d=[6 1 2 7 3 6 6 6 2 2 5]
%c=[1 2 1 2 1 3 5 4 3 4 8]
In=sparse(c,1,1);
Out=sparse(d,1,1);
nB=max(size(In,1),size(Out,1));
In(nB,1)=0; % force In to be of size nB-by-1
Out(nB,1)=0; % force Out to be of size nB-by-1
Out0=find(Out==0);  % vertices with 0 edges going out
d_=[d(:); Out0]; c_=[c(:); Out0]; % add necessary extra edges to permit the linking process
In(Out0)=In(Out0)+1; Out(Out0)=1; % updating In and Out with the extra added edges in d_, c_, so that now we have In=sparse(c_,1,1);Out=sparse(d_,1,1);
Incum=[0; cumsum(In(1:end-1))];
Outcum=[0; cumsum(Out(1:end-1))];
[cs,sc]=sort(c_);
[~,sd]=sort(d_);

%mod((1:numel(sd))'-Incum(cs)-1,Out(cs))+1 +Outcum(cs)
phi=zeros(size(sc));
phi(sc)=sd(mod((1:numel(sd))'-Incum(cs)-1,Out(cs))+1+Outcum(cs));

isequal(c_, d_(phi))


case '29May2024'
% the signature has been changed
% function [phi,theta,in,out]=linkdigraph(source,target,nE)
% [phi,theta,in,out]=linkdigraph(source,target,nE)
% Created 29May2024, last version 5May2025
% linkdigraph24a.m 29May2024 - can be improved with weighted vertices as input to be used with sort in the linking process
% to be used with the new version of linkdigraph24.m
% [s,t,u,Adj,Edges]=digraph24(source,target,nE)
%
% Given a graph sourse,target:E--V with E=(1:nE), calls linkdigraph24.m to make sure that a link (phi) is possible, and returns that link (phi) such that t=s(phi), with t and s obtained from 
% [s,t,u,Adj,Edges]=digraph24(source,target,nE)
% it also returns (theta) which is a cyclic order around each vertex, that is, a bijection with s(theta)=s

try
[s,t,u,Adj,Edges]=digraph(source,target,nE);
%[s,t,u,Adj,Edges]=digraph24(source,target,nE);
catch
[s,t,u,Adj,Edges]=digraph(source,target,size(target,1));
%[s,t,u,Adj,Edges]=digraph24(source,target,size(target,1));
end
nB=max(t); B=(1:nB);
nA=numel(t); A=(1:nA);
phi=zeros(nA,1);
theta=zeros(nA,1);
in=full(sparse(t,1,ones(size(t))));
out=full(sparse(s,1,ones(size(s))));
for b=1:nB
    Atb=A(t==b);
    Asb=A(s==b);
    if ~isempty(Asb)
        nRM=ceil(numel(Atb)/numel(Asb));
        RM=repmat(Asb,1,nRM);
        phi(Atb)=RM(1:numel(Atb));
        theta(Asb)=[Asb(2:end) Asb(1)];
    end
end %for
end % switch version
end % function
