function [phican,sphi]=linkcanonical(phi)
% [phican,sphi,sphinv]=linkcanonical(phi)
%
% 17Jun2025 - first trial
% 
% given an endomap phi, returns phican such that issorted([orb,ord]) with [orb,ord]=orbits(phi), it works well for invertible phi's

idom=@(x) (1:numel(x))'; % identity of domain of vector x considered as a morphism in Matlab-Octave category
[orb, ord, psi, deg, init, term, prin, conn] = orbits(phi);
SM=sortrows([orb ord idom(phi)])
% a more accurate version when phi is not an iso would be perhaps
%SM=sortrows([conn(orb) term(orb) init(orb) orb ord idom(phi)])
sphi=SM(:,end); % the permutation that makes phi sorted as desired
sphinv=zeros(size(phi)); % inicialize the inverse
sphinv(sphi)=idom(sphi); % compute the inverse
phican=sphinv(phi(sphi));

%{
% Examples
phi=[3 4 2 1]', pos=[0 1+1i 1 1i]'
[phican,sphi]=linkcanonical(phi)
linkdraw(pos(sphi),phican)
% disp([idom(phi) phi orb ord sphi phi(sphi) sphinv phican])

phi=randperm(25)
pos=homology(phi);
subplot(1,2,1)
linkdraw(pos,phi)
[phican,sphi]=linkcanonical(phi);
subplot(1,2,2)
linkdraw(pos(sphi),phican)

%}