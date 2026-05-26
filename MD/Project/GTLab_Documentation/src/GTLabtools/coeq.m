function [q,p1,p2]=coeq(d,c)
% [q,p1,p2]=coeq(d,c)
%
% given a graph (d,c), considered as two maps d,c:A-->B with A=(1:nA)' and B=(1:nB)', and such that d is surjective, returns q=coeq(d,c) the categorical coequalizer of d and c, so that q(d)==q(c) and for every  f:B-->C with f(d)==f(c) there exists a unique fbar:(1:numel(u)')-->C such that f=fbar(r(q)) for a unique fbar obtained as fbar=f(r) with [u,r,q]=unique(q); this replaces the output q indexed by representatives into a new q with linear indexes, so that the old q is recovered as u(q).
% the output values p1,p2 are such that if R=sparse(p1,,p2) then R+R' is the kernel pair relation of q 

nB=max(d); B=(1:nB)';
if ~isequal(unique(d),B)
    fprintf('Warning: d is not surjective, q is not guaranteed to be a coequalizer \n');
end

R0=sparse(d,c,1,nB,nB);
R1=triu(R0+R0');
R=R1; % Get the transitive closure of R1
Rtran=logical(R+R*R1);
while ~isequal(R,Rtran)
    R=Rtran;
    Rtran=logical(R+R*R1);
end

[s,t]=find(R+eye(nB));
E=[s,t]; Erows=sortrows(E);
p1=Erows(:,1); p2=Erows(:,2);
q(p1)=p2; q=q(:);

end % of function

% Examples
%{
d=[1 4 3 6 2]'; c=[3 5 2 5 1]';
[q, p1, p2]=coeq(d,c); 
[d c], [(1:max(d))' q], [(1:numel(p1))' p1 p2]


d=[1 4 3 6 2 5]'; c=[3 5 2 5 1 4]';
[q, p1, p2]=coeq(d,c); 
[d c], [(1:max(d))' q], [(1:numel(p1))' p1 p2]

d=[8 4 3 6 2 5  7 1]'; c=[7 5 2 5 1 4 8 1]';
[q, p1, p2]=coeq(d,c); 
[d c], [(1:max(d))' q], [(1:numel(p1))' p1 p2]

d_=(mod(1:5000,200)+1)'; c_=(mod([2 : 5000 1],200)+1)';
[q, p1, p2]=coeq(d,c); 
[d c], [(1:max(d))' q], [(1:numel(p1))' p1 p2]


%}