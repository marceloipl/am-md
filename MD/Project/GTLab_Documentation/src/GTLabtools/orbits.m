function [orb, ord, psi, deg, init, term, prin, conn] = orbits(phi)
% [idn, phi, orb, ord, psi, deg, ini, term, prin, conn] = orbits(phi)
%%
% Date: 2Fev2024
% % given a vector phi, interpreted as an endomap from the set of indexes 1:numel(phi)
% into itself, returns vectors of the same size as phi with the orbits (orb),
% the order inside each orbit (ord), a pseudo-inverse to phi (psi), and the degree
% of each point (deg); it also returns a vector with the initial points
% (init), the terminal points associated with each initial point, that is,
% the last point in oder of its principal component, and the level of the
% component (prin), 1 means principal, 2 is for secondary, 3 terciary etc
%
% (a) the orbit is a label which identifies the points a, phi(a), phi^2(a),...,phi^k(a),
% until a point phi^(k+1)(a) is found which is already a member of the current
% sequence or has already been considered in a previous sequence; it is initiated
% with the numbers 1, 2, 3,... up to the total number n0 of initial
% points (the ones with degree 0); the sequence continues with the numbers
% n0+1, n0+2, n0+3,...
% for the cyclic components which are characterized by the fact that each point
% in it is a regular point (degree -1)
%
% (b) the  order is the number k described in the previous item
%
% (c) pseudo-inverse to phi in the sense that it is consistent with the choice
% of principal components
%
% (d) degree in the sense that deg(a)=p means that there are p elements x such
% that phi(x)=a; the points in a cyclic component are labeled -1
%
% Example
% id=@(x) [1:numel(x)]';
% phi=[4 5 5 7 10 8 1 10 10 11 9 8 13 13];
% [ orb, ord, psi, deg, init, term, prin, conn] = orbits_wv04(phi);
% disp([id(phi) phi(:) orb ord psi deg ])
%      1     4     6     1     7    -1
%      2     5     1     1     2     0
%      3     5     2     1     3     0
%      4     7     6     2     1    -1
%      5    10     1     2     2     2
%      6     8     3     1     6     0
%      7     1     6     3     4    -1
%      8    10     3     2     6     2
%      9    10     1     5    11     1
%     10    11     1     3     5     3
%     11     9     1     4    10     1
%     12     8     4     1    12     0
%     13    13     5     2    14     2
%     14    13     5     1    14     0
%
% disp([init; term; prin ; conn ])
%      2     3     6    12    14
%      9     3     8    12    13
%      1     2     2     3     1
%      1     1     1     1     5
%
% sortrows([orb ord id(orb)])
%      1     1     2
%      1     2     5
%      1     3    10
%      1     4    11
%      1     5     9
%      2     1     3
%      3     1     6
%      3     2     8
%      4     1    12
%      5     1    14
%      5     2    13
%      6     1     1
%      6     2     4
%      6     3     7

% given a vector phi
try
    phi=phi(:); % force phi to be a column vector
    phi(phi); % make sure that phi can be interpreted as an endomap
catch ME
    disp('The input vector phi must be interpreted as an endomap from the set of indexes 1:numel(phi) to itself:')
    error(ME.message)
end
n=numel(phi);
idn=(1:n)'; % the identity vector of the same size as phi
% initializing outputs
orb=zeros(n,1);
ord=zeros(n,1);
psi(phi)=idn; psi=psi(:); % ensure column vector
deg=repmat(-1,n,1);  % the degree is initiated by default as cyclic
% determining initial points (those with degree 0)
init=idn(~ismember(idn,phi))'; % vector with indexes of initial points
m=numel(init); % number of initial elements
term=zeros(1,m); % the vector term will contain the terminal point of each initial point with respect to the choice of the principal components
prin=zeros(1,m); % the vector prin will contain the principal components, secondary, terciary, etc
conn=zeros(1,m); % the connected component, relative to the principal component

% the iterative process
for j=1:numel(init)
    a=init(j); k=1;
    orb(a)=j;
    ord(a)=k;
    psi(a)=a; % the initial points become fixed points when inverted
    deg(a)=0;
    while ~orb(phi(a))   % phi(a) is not yet a junction point
        psi(phi(a))=a; % is defined before a is uptaded as a=phi(a)
        a=phi(a); k=k+1;  % a is updated to phi(a) and k to k+1
        orb(a)=j;
        ord(a)=k;
        deg(a)=1;
    end
    % phi(a) is a junction point (orb(phi(a)) is not zero)
    if orb(a)==orb(phi(a))  %  principal component
        deg(phi(a))=deg(phi(a))+1;
        term(j)=a;
        prin(j)=1;
        conn(j)=j;  % the same as orb(a)
    else % secondary component
        deg(phi(a))=deg(phi(a))+1;
        term(j)=a; % each secondary component will use the terminal point of its principal component
        prin(j)=prin(orb(phi(a)))+1;
        conn(j)=conn(orb(phi(a)));
    end
end

% It remains to compute the cyclic orbits
% e.g. phi=[2 3 4 5 1 1 6 6 1 3 12 13 14 15 12 13 16 17]'; phi=[phi' 20 21 22 19 23 25 24]; disp([idphi phi orb ord psi])
q=m;  % m=numel(init), initialize component label, in oder to continue the non-cyclic components
while any(~orb)
    a=find(~orb,1);
    q=q+1;  % update component lablel
    orb(a)=q;  % set component label
    ord(a)=1;  % initialize counting the order of current component
    %deg(a)=-1;  % we have choosen to single out the cyclic points by pre-setting deg=-1
    while ~orb(phi(a))
        a=phi(a);
        orb(a)=q;
        ord(a)=ord(psi(a))+1;
        deg(a)=-1;
    end
end

