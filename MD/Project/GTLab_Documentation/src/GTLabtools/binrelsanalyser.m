function [prop,Rtran,Ranti]=binrelsanalyser(R)
% [prop,Rtran,Ranti]=binrelsanalyser(R)
%
% 30May: companion to endomaps_exe.m
% given a sparse adjacency matrix R of a binary relation, ckecks for the properties: reflexivity, transitivity, symmetry, anti-symmetry, and returns the new matrix Rtran with the transitive closure if the relation R happens not to be transitive as well as a new relation Ranti if the relation R (in case it is not an equivalence relation) is not antisymmetric. If the relation R is transitive then Rtran is empty and if R is symmetric then Ranti is also empty. Only in the case when R is neither symmetric nor antisymmetric does the relation Ranti contain the result of forcing antisymmetry by first defining the equivalence relation x~y iff R(x,y) and R(y,x), so that Ranti is, in some sense, R/~ which is anti-symmetric with adjacency Ranti
% prop is a number from 0 to 5 with the following interpretarion:
% prop=0, R is not a square matrix
% prop=1, R is reflexive but not transitive
% prop=2, R is a preorder (reflexive and transitive), but not and equivalence relation (not a partial order)
% prop=3, R is an equivalence relation
% prop=4, R is a partial order (reflexive, transitive and anti-symmetric)
% prop 5, R is not an equivalence relation but is a preorder (reflexive and transitive) and Ranti is the biggest partial order contained in R
% prop 6, otherwise
%
% testing examples at the end of the file

prop=6; Rtran=[]; Ranti=[];
[n1,n2]=size(R); idV=(1:n1)'; % the set of indexes for the vertices V={1,...,n1} if R \subseteq V \times V is seen as a graph with s,t:nonzeros(R)-->V defined by [s,t]=find(R)

if ~isequal(n1,n2),
    disp('R is not a square matrix')
    prop=0;
    return
end
if ~isequal(diag(R),ones(n1,1))
    disp('R is not relexive')
    prop=0;
else
    disp('R is reflexive')
    prop=1;
end

if ~(all(all(logical(R^2)<=R))) && prop==1
    disp('R is not transitive, the output Rtran is its transitive closure')
    Rn=R; % Get the transitive closure of R (see also Warshall Algorithm)
Rtran=logical(Rn+Rn*R);
    while ~isequal(Rn,Rtran)
        Rn=Rtran;
        Rtran=logical(Rn+Rn*R);
    end    
elseif prop==1
    prop=2;
    disp('R is reflexive and transitive (also called a preorder)')
end

if isequal(R,R')
    disp('R is symmetric')
    if prop==2, prop=3; disp('R is an equivalence relation') end
elseif ~isequal((R & R'),eye(size(R)))
    disp('R is not symmetric (hence not an equivalence relation) nor anti-symmetric (hence not a partial order)')
    if prop==2 
        prop=5;
        disp('but it is a preorder (reflexive and transitive), so, forcing anti-symmetry by identifying x~y iff R(x,y) and R(y,x), a partial order is obtained,')
    disp('Ranti is the adjacency matrix of the partial order R/~')
    [s,t]=find(R & R');   % identifying those points (x,y) for which anti-symmetry fails, by setting x~y iff R(x,y) and R(y,x)
q=coequalizer(s,t);  % the quotient map for the induced partition from the equiv rel (R & R')
Ranti=R(q==idV,q==idV);
end
else
    disp('R is anti-symmetric')
    if prop==2, prop=4; disp('R is a partial order'), end
end


%{
Testing examples
k=3, n=2*k, d=1:n; c=mod(k*d,n)+1; R=sparse(d,c,1,n,n)
[prop,Rtran,Ranti]=binrelsanalyser(R)

% not square R
R=[0 1 0; 0 0 1]; [prop,Rtran,Ranti]=binrelsanalyser(R)

% not reflexive but square
R=[0 1 0; 0 0 1; 0 0 0]; [prop,Rtran,Ranti]=binrelsanalyser(R)

% reflexive but not transitive
R=[1 1 0; 0 1 1; 0 0 1]; [prop,Rtran,Ranti]=binrelsanalyser(R)

% reflexive and transitive (or preorder)
R=[1 1 1; 1 1 1; 0 0 1]; [prop,Rtran,Ranti]=binrelsanalyser(R)

% symmetric but nothing else
R=[0 1 0; 1 0 0; 0 0 1]; [prop,Rtran,Ranti]=binrelsanalyser(R)

% an equivalence relation
R=[1 1 0; 1 1 0; 0 0 1]; [prop,Rtran,Ranti]=binrelsanalyser(R)

% a partial order
R=[1 1 0; 0 1 0; 0 0 1]; [prop,Rtran,Ranti]=binrelsanalyser(R)

%}