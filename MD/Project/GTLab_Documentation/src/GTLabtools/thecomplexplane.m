function [F,i2s,s2i,nA]=thecomplexplane(nB,nX)
% [F,i2s,s2i,nA]=thecomplexplane(nB,nX)
% 16Mar2025
%
% Given two numbers nB and nX, consider B=1:nB, X=1:nX and A=1:nB^nX with the bijections i2s:A-->B^X and s2i:B^X-->A, with i an index in 1:nA and s a vector of subindexes representing a map from X to B; construct the function F:A-->ComplexPlane as in dimagmaenergy.m

nA=nB^nX;
% the bijection
i2s=@(i) mod(floor(bsxfun(@rdivide,i(:)-1,nB.^((1:nX)-1))),nB)+1;
s2i=@(s) (s-1)*(nB.^((1:nX)'-1))+1;

% the auxiliar maps to define F

% construct the map f:X-->Integers
f=zeros(1,nX); % inicialize f with zeros
f(2:2:end)=-(2:2:nX)/2; % f(x) if x is even
f(1:2:end)=((1:2:nX)-1)/2 % f(x) if x is odd
% construct the map g:B-->Complex
g=[exp(2*pi*1i*(nB-(1:nB-1))/(nB-1)), 0]

% construct the map F:A-->Complex
F=@(s) sum(g(s).*3.^repmat(f,size(s,1),1),2);

% Examples of application

[F,i2s,s2i,nA]=thecomplexplane(7,3);
plot(F(i2s(1:nA)),'.')

