function [g,Hc,Hnc, sc, snc]=homology(phi,n)
% [g,Hc,Hnc]=homology(phi,n)
%
% 16Jun2025, added n as input for returning Hnc with the number of points up to degree n; that is, Hnc(i,3+n)= number of points on the non-cyclic connected component i with degree n [for the new output sc and snc see lines 31 and 75]
% 24May2024 - simplified and tested (the old version is homology_.m)
% given an endomap phi:X-->X with X=(1:numel(phi)) computes the orbits of phi and returns a complex vector g of the same size as phi with geometic realization together with Hc an homology characteristic vector with the number of elements is each cyclic component as well as Hnc a homology matrix with information on each non-cyclic connected component
%
% is an alternative to orbits24.m, orbits24_inprogress and to homologymatrix.m (and to homology_.m)
%
% Example (see also at the end of file)
% phi=[2 5 5, 15 15 16, 9 9 10, 11 8 13, 14 12 19, 21 16 19, 19 21 17, 17 ]
% [g,Hc,Hnc]=homology(phi)
% linkdraw23(phi,g,1:numel(phi))
%
% history: last stabe version24May; bug fixed on line 33 (31May2024)

[orb,ord,psi,deg,init,term,prin,conn]=orbits(phi);
idA=(1:numel(phi))';
g=zeros(size(phi)); % initialize g, the geometrical realization of phi

% the case of cyclic components
iscyc=(deg==-1);
cyc1=find(iscyc & ord==1); % indexes for first cycles
Hc=sortrows([ord(psi(cyc1)), orb(cyc1)]); % the second row contains the reference to the orbit, should be discarded for homology comparisons
% Example: >> phi=[2 3 4 1 6 7 5 9 8]; 
% Hc =[ 2   3   4]'

% the geometric realization for the cyclic components
SM=sortrows([orb(iscyc), ord(iscyc), idA(iscyc)]);
g(SM(:,3))=(1i*2*pi*linspace(1/nnz(iscyc),1,nnz(iscyc)));
sc=SM(:,3); % permutation that turns phi into a canonical phi by phican=phi(sc), added as output subrepticially


% the case of non-cyclic components
connop=conn'; prinop=prin'; % fixing a bug 31May2024
SM=sortrows([connop(orb(~iscyc)), prinop(orb(~iscyc)), ord(phi(term(orb(~iscyc)))), orb(~iscyc), ord(~iscyc), idA(~iscyc)]);  % sorting the data to be used further on
[~,r,q]=unique(conn); 
% conn is replaced by q with linear indexes and r(q)=conn'
Prin=sparse(q,prin,ord(term)); % combines each index of connected components with each one of its principal secondary tertiary components etc

% version from 24May2024
% Hnc=zeros(numel(r),6); Hnc(:,6)=r; % the 6th column is only for reference
% the total number of points in each unique connected component
% Hnc(:,1)=full(sum(Prin,2));
% the period of each unique connected component
%Hnc(:,2)=ord(term(r))-ord(phi(term(r)))+1;
% the degrees of each connected component
%for i=1:numel(r)
%isconni=(conn(orb(~iscyc))==r(i))';
%Hnc(i,3)=nnz(isconni & deg(~iscyc)==0); % initial points
%Hnc(i,4)=nnz(isconni & deg(~iscyc)==2); % junction points with deg==2
%Hnc(i,5)=nnz(isconni & deg(~iscyc)>=3); %junction points with deg>=3
%end
% New version 16Jun2025
try n, catch n=0, end

Hnc=zeros(numel(r),3+n); Hnc(:,end)=r; % the last column is only for reference
% the total number of points in each unique connected component
Hnc(:,1)=full(sum(Prin,2));
% the period of each unique connected component
Hnc(:,2)=ord(term(r))-ord(phi(term(r)))+1;

for i=1:numel(r)
isconni=(conn(orb(~iscyc))==r(i))';
for j=0:n % number of homology that counts the degree of each point in a non-cyclic component i
Hnc(i,3+j)=nnz(isconni(:) & deg(~iscyc)==j); % initial points n=0
%Hnc(i,3+1)=nnz(isconni & deg(~iscyc)==1); % connecting points with deg==1
%Hnc(i,3+2)=nnz(isconni & deg(~iscyc)==1); %junction points with deg==2
%Hnc(i,3+3)=nnz(isconni & deg(~iscyc)>=3); %junction points with deg=3
%end
% the total number of points with degree >= n is obtained by summing all degress taking the difference to the total number of points in the component wich is stored at Hnc(i,1)
end
end

snc=SM(:,3); % permutation that turns phi into a canonical phi by phican=phi(sphi), added as output subrepticially

% from here on is the version from 24May2024

% computing the starting and ending angles for each unique connected component
th=2*pi*(cumsum([0; Hnc(:,1)]))/sum(Hnc(:,1)); 
% each component starts at th(q(i)) and ends and th(q(i)+1) angles between 0 and 2pi

% computing the geometry of each unique connected component associated with its principal, secondary etc, sub-components
[I,J]=find(Prin); 
% each pair (I(i),J(i)) is an index of a connected component (I(i)) combined with an index of a principal, secondary etc sub-component J(i); the component indexed by I(i) is r(I(i)) 
for i=1:numel(I)
    iscompi=(SM(:,1)==r(I(i)) & SM(:,2)==J(i));
        ni=nnz(iscompi);
        thi=th(I(i));
        rangi= 0.5 *(1 + linspace(-1+1/ni,1-1/ni,ni));
        rangethi=min(pi/2,th(I(i)+1)-thi);
        g(SM(iscompi,6))=log(sqrt(J(i)+1)) + 1i*(thi+ rangi*rangethi);
end

g=exp(g);

%Hnc=sortrows(Hnc(:,1:5)); % sorting and discarding the 6th column
Hnc=sortrows(Hnc(:,1:end-1)); % sorting and discarding the last column
Hc=Hc(:,1); % removing the 2nd column in Hc which was just a reference to the cyclic component 

end % function
%% More examples
%{
disp('Example 1')
phi=[2 5 5, 15 15 16, 9 9 10, 11 8 13, 14 12 19, 21 16 19, 19 21 17, 17 ]
[g,Hc,Hnc]=homology(phi)
linkdraw23(phi,g,1:numel(phi))

disp('Example 2')
phi=[2:20, 5, 22:30, 25, 32:45, 41, 47 48 46];
[g,Hc,Hnc]=homology(phi)
linkdraw23(phi,g,1:numel(phi))

disp('Example 3')
phi=[2 2 4 5 6 2 8 9 2]
[g,Hc,Hnc]=homology(phi)
linkdraw23(phi,g,1:numel(phi))

disp('Example 4')
phi=[2 3 4 1 6 7 5 9 8];
[g,Hc,Hnc]=homology(phi)
linkdraw23(phi,g,1:numel(phi))

% this one gives an error
homology([2 2 1])

%added 16Jun2025
phi=[2 5 5, 15 15 16, 9 9 10, 11 8 13, 14 12 19, 21 16 19, 19 21 17, 17 ]; 
phi=[phi numel(phi)+[2:20, 5, 22:30, 25, 32:45, 41, 47 48 46]];
phi=[phi,phi]; 
phi=[phi,phi];
phi=[phi,phi];
[g,Hc,Hnc]=homology(phi,8)
linkdraw(g,phi)

% still problems with
[g,Hc,Hnc]=homology([2 2 1],4)
linkdraw(g,[2 2 1])

% to get a canonical phican
phi=10:-1:1
phi=[3 1 2 5 6 4]
phi=[7 1 3 6 2 4 8 5]
phi=randperm(20)
[g,Hc,Hnc, sc, snc]=homology(phi)
linkdraw(g,phi)
%phican=phi([snc;sc]) %not quite there yet, fixed as follows
sphi=[snc;sc]
sphinv(sphi)=(1:numel(sphi))'
phican=sphinv(phi(sphi))

%}