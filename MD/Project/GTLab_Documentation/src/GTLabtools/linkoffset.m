function gh=linkoffset(h,g,phi,nA)
% Given a link structure A=1:nA, phi, g with only one cyclic component and an offset value h, returns the new coordinates gh of the link structure. If h is negative the offset goes into the interior of the region defined by the curve, if h is positive then the offset goes to its exterior region. The connectivity of the link, that is, phi, does not change
% Further details in Caderno 146 - 10Feb2025

try
idA=(1:nA)'; 
catch
idA=(1:numel(g))';
end
phi=phi(:);

% Swicht implementation; New=10Feb2025; Old=before 2024
%Implement='Old';
Implement='New';


switch Implement
case 'New'  % 10Feb2025
pos=g(:);
dir=pos(phi)-pos;
alpha=angle(dir(phi)./dir);
rho=abs(dir);

rho_h=zeros(size(phi));
dir_h=zeros(size(phi));
pos_h=zeros(size(phi));

rho_h(phi)=h*(tan(alpha/2)+tan(alpha(phi)/2))+rho(phi);

ace_h=log(rho_h(phi)./rho_h)+1i*alpha;
dir_h(phi)=exp(cumsum(ace_h));
pos_h(phi)=cumsum(dir_h);
%clf, linkdraw(pos); hold on, linkdraw(pos_h); hold off

q=coeq(idA,phi);
%gh=pos_h + h*(1+1i*tan(alpha(q)/2));
gh=pos_h;

case 'Old' % Before 2024 ---------------------
A=idA;
d=g(phi)-g; dA=d(A);
dd=d(phi)./d; ddA=dd(A);
phiA=phi(A); gA=g(A);

alpha=angle(ddA);
dh=abs(dA); % initialize dh, where dh means dA after being offseted by h

%h=-1.8;
dh=dh+h*tan(alpha/2);
dh(phiA)=dh(phiA)+h*tan(alpha/2);

dh=dh.*exp(sqrt(-1)*angle(dA)); % the angles do not change with offset

offset=h/cos(alpha(end)/2)*exp(sqrt(-1)*((pi + angle(dA(1)) - alpha(end)/2)+pi/2));
gh=cumsum([gA(1)+offset, dh ]);
disp(abs(g(1)-g(end)))
gh(end)=[];
end

% Examples:
 %{
 g=[1 1i -1 -1i]; phi=[2 3 4 1];
gh=linkoffset(.5,g,phi);
linkdraw(g,phi); hold on, linkdraw(gh,phi);
 %}
 
 %{
 connsig=[3 3 3 3]  % connected comp signature
 n=[0 cumsum(connsig)];
conn=cumsum(ismember((1:n(end))',n(1:end-1)+1));
phi=[2:n(end) 0]'; % auxiliar step
phi([n(2:end)])=n(1:end-1)+1;
idA=(1:n(end))';
q=coeq(idA,phi);
orb=orbits(phi);
linkdraw(homology(phi),phi)
g=homology(phi);
gh=linkoffset(.5,g,phi)
 %}