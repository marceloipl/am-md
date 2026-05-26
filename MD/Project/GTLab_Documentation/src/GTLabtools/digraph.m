function [s,t,u,Adj,Edges]=digraph(source,target,nE)
% [s,t,u,Ajd]=digraph(sourse,target)
%
% given a graph structure (source,target,nE) in which source and target are matrices with the same number of columns and at least nE lines interpreted as a directed graph with nE edges indexed in 1:nE and such that an edge i in 1:nE starts at vertex source(i,:) and ends at vertex target(i,:) with the lines in source and target interpreted as vectors in euclidean space
% Returns the indexed digraph (s,t,u) with numel(s)<=nE the number of edges i in i:nE with the property that there exists at least one j in 1:nE such that source(j,:)==target(i,:), the matrix u contains the unique rows of source and Adj is the adjacency matrix of the digraph (s,t,nA) with nA=numel(s)=numel(t) such that Adj(i,j)=k has the meaning that k in 1:nE is the original index of edge k in 1:nE in the sense that source(k)=s(i) and target(k)=t(j) 
% Edges is the adjacency matrix of (source,1:nE)  where Edges(i,j) can be 2, 1,or -1  accordingly to: Edges(i,j)=1, the original edge j with domain i has the property described before ; -1 does not have the property above; 2 is a chosen edge for vertex index b in B=1:size(u,1) 
% Example see also the end of file
% source=[1 2; 2 3; 3 4 ]
% target=[3 4; 3 4; 1 2; 3 2]
% [s,t,u,Adj,Edges]=digraph24(source,target,3)

if ~isequal(size(source,2),size(target,2)), error('LinksToolbox error: source and target must have the same number of columns'), end
try
E=(1:nE)'; 
dom=source(E,:);
cod=target(E,:); disp('yes')
catch 
E=(1:min(size(source,1),size(target,1)))'; 
dom=source(E,:);
cod=target(E,:);
fprintf('Linkstoolbox: the number of edges considered was: %d \n ', numel(E))
end
    
[u,e_,d]=unique(dom,'rows');
[k,c]=ismember(cod,u,'rows');
s=d(k);
t=c(k);
Adj=sparse(s,t,E(k));
cval=zeros(size(k)); cval(k)=1; cval(~k)=-1; cval(e_)=2;
Edges=sparse(d,E,cval);

end  % end of function
disp(example)
source=[1 2; 2 3; 3 4; 2 3 ; 3 4]
target=[3 4; 3 4; 1 2; 3 2 ; 2 3]
[s t u Adj,Edges]=digraph24(source,target)

% generating a cube
u=fillzeros([0 0 0]',0:1)';
dc=fillzeros([0 0]',1:8)'; d=dc(:,1); c=dc(:,2);
uu=u(c,:)-u(d,:);
k=dot(uu',uu')==1;
s=d(k); t=c(k);
source=u(s,:); target=u(t,:);
[s_ t_ u_ Adj,Edges]=digraph24(source,target);
disp([s t s_ t_])
disp([u, u_])
full(Adj)