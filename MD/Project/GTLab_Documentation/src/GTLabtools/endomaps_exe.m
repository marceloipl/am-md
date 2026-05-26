function [R,Ranti]=endomaps_exe(rel,nX)
% R=endomaps_exe(rel,nX)
%
% 30May: this is an executable version of endomaps in which nX is provided as well as the number of the relation (rel) to perform the calculations for R1 to R6 (so that rel must be a number from 1 to 6)
% the default values are rel=1 and nX=3

try
nX < 9; rel < 7
catch
nX=3; rel=1; warning('Input nX was set to 3 and input (rel) was set to 1')
end

% initializations, for more details see endomaps.m of the documentation ClassroomMDPL24.pdf section 8

%nX=3;  % it is now prepared to be scaled with n = 4,5,...
eyeX=(1:nX);  % the identity endomap on X={1,...,nX}
nB=nX^nX;
B=(1:nB)';  % the set of all endomaps EndX linearly indexed from 1 to n^n

% the composition of endomaps as if they were row vectors
composerows=@(g,f) g(sub2ind(size(f),repmat((1:size(f,1))',1,size(f,2)),f));

% the bijection between size nB and EndX
i2f=@(i) mod(floor(bsxfun(@rdivide,i(:)-1,nX.^((1:nX)-1))),nX)+1; % this is a more efficient alternative for greater nX
%EndXop=fillzeros(zeros(size(eyeX)),eyeX)';  % note that we are now working with rows instead of columns
%i2f=@(i) EndXop(i,:);  % obsolete
f2i=@(f) (f-1)*(nX.^((1:nX)'-1))+1;

% ckecking that i2f and f2i are inverse to each other, without using EndXop, here and elsewhere EndXop is being replaced by i2f(B)
isequal(f2i(i2f(B)),B)
isequal(i2f(B),i2f(f2i(i2f(B))))

% finding invertible and monotone endomaps
[fsorted,finvert]=sort(i2f(B),2);
phi=find(f2i(fsorted)==B); % indexes to monotone endomaps
alpha=find(all(fsorted==repmat(eyeX,nB,1),2)); % indexes to invertible endomaps
alphainv=f2i(finvert(alpha,:)); % indexes to the respective inverses of alpha

% ckecking that indexes in alpha and alphainv correspond to each other
isequal(composerows(i2f(alpha),i2f(alphainv)),repmat(eyeX,numel(alpha),1))
isequal(composerows(i2f(alphainv),i2f(alpha)),repmat(eyeX,numel(alpha),1))

% computing the relation for the given case (rel)
R=binaryrel(rel,nX,phi,alpha,alphainv);

% post processing the relation
% is it an equivalence relation?
% if not an equivalence relation then construct the Hasse diagram of the larger partial order contained in R
% this analysis if performed by the function binrelsanalyser.m
[prop,Rtran,Ranti]=binrelsanalyser(R);


if prop==3
    % is an equivalence relation
fprintf('The relation R%d is an equivalence relation\n',rel)    
[s,t]=find(R);
qR=(coequalizer(s,t));
disp('Number of equivalence classes')
numel(unique(qR))   
%SM=sortrows([qR B i2f(B)]) % i2f(B) is  EndXop
%showmethemaps(SM)
%    fi=logical(sum(R-eye(size(R)),2))
%    [s t]=find(R-diag(fi));
%[phi,theta,in,out]=linkdigraph24a(s,t);
g=homology(qR);
linkdraw24(g,qR,(1:numel(qR))',0,0);

elseif prop==5
    % Ranti is a partial order, let us try to visualize its Hasse diagram
    fprintf('In correcting the failure of anti-symmetry for the relation R%d the following partition is created\n', rel)
    [s,t]=find(R & R');
qR=(coequalizer(s,t));
disp('Number of equivalence classes')
numel(unique(qR))   
g=homology(qR);
subplot(2,1,1)
linkdraw24(g,qR,(1:numel(qR))',0,0);
%pause
title('correction for the failure of anti-symmetry')
    % [s,t]=find(Ranti);  
    % [s,t]=find(Ranti-eye(size(Ranti)));
    fi=logical(sum(Ranti-eye(size(Ranti)),2));
    [s t]=find(Ranti-diag(fi));
[phi,theta,in,out]=linkdigraph24a(s,t);
g=homology(phi);
subplot(2,1,2)
linkdraw24(g,phi,s,0);
%R=Ranti; % return as output the antisymmetric relation
end






function R=binaryrel(rel,nX,phi,alpha,alphainv)
% auxiliar function to compute the adjacency matrix of each one of the 6 relations

% the preliminaries
nB=nX^nX;
B=(1:nB)';
composerows=@(g,f) g(sub2ind(size(f),repmat((1:size(f,1))',1,size(f,2)),f));
i2f=@(i) mod(floor(bsxfun(@rdivide,i(:)-1,nX.^((1:nX)-1))),nX)+1; 
f2i=@(f) (f-1)*(nX.^((1:nX)'-1))+1;

% the 6 cases
switch rel
case 1 % the first relation 
% The first relation: R=Til
I=alpha; J=alphainv; % B=(1:nB) as defined above
I_=(1:numel(I))'; % indexes to I
[i,b]=ndgrid(I_,B); 
j=J(i);   % this ensures that J depends on I 
i=I(i);   % this replaces i from 1:6 by the values of I=alpha
Fi=i2f(i); % convert indexes into maps as rows in EndXop
Fb=i2f(b);
Fj=i2f(j);  
m=f2i(composerows(Fi,composerows(Fb,Fj))); % compute Fi(Fb(Fj)) and convert the result to an index
% the desired relation is:
R=logical(sparse(b(:),m(:),1));
%------------------------------------------------------
case 2
% The second relation: R=Til_L
J=alpha; % B=(1:nB) 
[b,j]=ndgrid(B,J);  % j is ranging all values in alpha 
Fb=i2f(b); % convert indexes into maps as rows in EndXop
Fj=i2f(j);  
m=f2i(composerows(Fb,Fj));
R=logical(sparse(b(:),m(:),1));
%------------------------------------------------------
case 3
% The third relation: R=Til_R
I=alpha; % B=(1:nB) 
[i,b]=ndgrid(I,B);  % i is ranging all values in alpha 
Fi=i2f(i); % convert indexes into maps as rows in EndXop
Fb=i2f(b);  
m=f2i(composerows(Fi,Fb));
R=logical(sparse(b(:),m(:),1));
%------------------------------------------------------
case 4
% The fourth relation: R4=LeqL
%[alpha,alphainv]=find(circ==f2i(Xop)); % invertible elements
%phi=B(all(EndXop==sort(EndXop,2),2));  % monotone maps
J=phi; % B=(1:nB) % phi is indexing the monotone maps
[b,j]=ndgrid(B,J);  % j is ranging all values in phi 
Fj=i2f(j); % convert indexes into maps as rows in EndXop
Fb=i2f(b);  
m=f2i(composerows(Fb,Fj));
R=logical(sparse(b(:),m(:),1));
%-------------------------------------------------------
case 5
% The fifth relation: R5=LeqR
%[alpha,alphainv]=find(circ==f2i(Xop)); % invertible elements
%phi=B(all(EndXop==sort(EndXop,2),2));  % monotone maps
J=phi; % B=(1:nB) % phi is indexing the monotone maps
[b,j]=ndgrid(B,J);  % j is ranging all values in phi 
Fj=i2f(j); % convert indexes into maps as rows in EndXop
Fb=i2f(b);  
m=f2i(composerows(Fj,Fb));
R=logical(sparse(b(:),m(:),1));
%-------------------------------------------------------
case 6
% The sixth relation: R6=Leq
%[alpha,alphainv]=find(circ==f2i(Xop)); % invertible elements
%phi=B(all(EndXop==sort(EndXop,2),2));  % monotone maps
I=phi; J=phi; % B=(1:nB) % phi is indexing the monotone maps
[i,b,j]=ndgrid(I,B,J);  % j is ranging all values in phi 
Fj=i2f(j); % convert indexes into maps as rows in EndXop
Fi=i2f(i);
Fb=i2f(b);  
m=f2i(composerows(Fi,composerows(Fb,Fj)));
R=logical(sparse(b(:),m(:),1));
end % switch
end % subfunction
end % main function
