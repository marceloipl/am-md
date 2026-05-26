disp("Aula MDPL25 dia 7Mar2025");

% implmenting fillzeros
A0=[0 2 0 -1]
A1=A0(:)
z=find(A1==0)
nz=numel(z)
[a b]=ndgrid([1 3], 1:3)
ab=[a(:) b(:)]
A=repmat(A1,1,size(ab,1))
A(z,:)=ab'
%{
% has a problem
for i=1:nz
    A(A1==-i,:)=A(z(i),:);
end
%}
A(A1<0,:)=A(z(abs(A1(A1<0))),:);
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
% replacing negative values by the respective counterparts
=======
=======
=======
<<<<<<< HEAD
% replacing negative values by the respective counterparts
=======
>>>>>>> origin/master
>>>>>>> origin/master
=======
>>>>>>> origin/master
% replacing negative values by the respective counterparts

% A is the collection of all endomaps of X=1:3
A=fillzeros([0 0 0]',1:3)
Circ=zeros(27,27);
for i=1:27
for j=1:27
ui=A(:,i);
uj=A(:,j);
uij=ui(uj);
Circ(i,j)=find(ismember(A',uij(:)','rows'));
end
end

% identify the neutral element:
% (a) in the index matrix A; (b) in the multiplication table Circ
find(ismember(A',[1 2 3],'rows'))
find(ismember(Circ,1:27,'rows'))

% identify the invertible elements (i.e. permutations):
% (a) in the index matrix A; (b) in the multiplication table Circ
ismember(sort(A)',[1 2 3],'rows')'
ismember(sort(Circ,2),1:27,'rows')'
% disp(A)

% another possibility is
Perm=perms(1:3)'
Inv=find(ismember(A',Perm','rows'))

% Define the relation of equivalence up to conjugation, i.e., for all u and v in A, we say that u is related to v (and write uRv) iff there exists p in Per such that p(u)=v(p)

R=zeros(27,27);
for i=1:27
for j=1:27
for k=1:6
u=A(:,i);
v=A(:,j);
p=Perm(:,k);
if p(u)==v(p), R(i,j)=1; end
end
end
end

% to extract the graph associated with the relation R we do:
[d,c]=find(R);

% And the quotient A/R is computed with the coequalizer
q=coeq(d,c)

% The number of equivalence classes is:
numel(unique(q))
<<<<<<< HEAD
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
>>>>>>> origin/master
>>>>>>> origin/master
>>>>>>> origin/master
=======
>>>>>>> origin/master
