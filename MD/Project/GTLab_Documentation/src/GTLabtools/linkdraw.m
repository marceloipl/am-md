function [hq,ht]=linkdraw(pos,phi,typeofdraw,places,labels,idA)
% [hq,ht]=linkdraw(pos,phi,typeofdraw,places,labels,idA)
%
% Dez2024 (last stable version, has several differences when compared with its predecessors, namely typeofdraw)
% this is an alternative version to linkdrawplusplus.m
% Given:
% 1- pos, a vector with positions of the vertices (either complex or 3d coordinates)
% 2 - phi, the linking map, so that a digraph (dom,cod) is recovered as dom=pos cod=pos(phi,:)
% 3- typeofdraw = 0 1 2 3 where:
        % 0 - no arrowhead, no labels
        % 1 -  arrowhead,  labels
        % 2 - arrowhead, no labels
        % 3 - no arrowhead, labels
% 4 - places, labels placement, a value between 0 and 1 indicating the placement of edge a in idA placed at position pos(a)+places*dir(a), where dir=pos(phi)-pos
% 5 labels, the values (associated to the indexes in idA) to be converted into string text that will be placed at the positions indicated by places 
% 6 - idA, a subset of indexes of edges indicating which ones are labeled, must be the same size as labels
% Remarks:
% If input (5) labels is not provided then num2str(idA) is used insatead with idA=(1:numel(phi))';
% If input (4) places is not provided then 0.5 is used
% If input (6) is not provided then (1:numel(phi)) is used
% If input (3) is not provided then 0 is used (arrowhead, nolabels)
% If input (2) is not provided then phi=[2:numel(pos) 1] is assumed to be cyclic
% the output hq is the handle to the quiver plot whereas ht is the handle to the text labels
%

%{ 
%Simple example:
 N=6; idA=(1:N)';
 phi=[2:N 1]'; pos=exp(idA/N*2*pi*1i);
 lindraw(pos,phi), lindraw(pos,phi,1)
 lindraw(pos,phi,0), lindraw(pos,phi,2), lindraw(pos,phi,3)
 lindraw(pos,phi,3,0)
 linkdraw(pos,phi,3,0,idA/6*360)
 linkdraw(pos,phi,3,0,idA/6*360,1:2:5)
%}

N=nargin(); % number of input arguments, to determine the behaviour of the function

hq=0; ht=0; drawdim=2; drawlabels=1; drawarrowheads=1;

% taking care of input (1) pos
try
[nlin, ncol]=size(pos);
if nlin==1 % assume it is a row vector of complex entries and convert it into a column vector
pos=pos(:);
nlin=ncol;
%drawdim=2; % drawing dimension 2d use quiver
elseif ncol==2 % assume it is a nlin-by-2 matrix with [x,y] entries by row and convert it into a complex column vector
pos=pos(:,1)+1i*pos(:,2);
%drawdim=2; % drawing dimension 2d use quiver
elseif ncol>=3 % assume it is nlin-by-3 matrix with [x y z] coordinates by row and use quiver3
drawdim=3; % drawing dimension 3d use quiver3
end
catch
error('Linkstoolbox: the input pos must be a complex vector or else a nlin-by-3 matrix')
end

% adjusting remaining inputs
switch N % number of inputs by user
case 1 % only pos is given
phi=[2:nlin,1]';
typeofdraw=1; % yes yes
places=.5;
idA=(1:nlin)';
labels=idA;
case 2 % pos and phi are given
typeofdraw=1; % yes yes
places=.5;
idA=(1:nlin)';
labels=idA;
case 3 % the type of draw if given
places=.5;
idA=(1:nlin)';
labels=idA;
case 4  % places are given
idA=(1:nlin)';
labels=idA;
case 5  % labels are  given
idA=(1:nlin)';
case 6 % all six inputs are given
labels=labels(idA);
labels=labels(:);
end


% defining dir the direction of each edge
dir=pos(phi,:)-pos(:,:);

place=pos(idA,:)+places.*dir(idA,:);

switch drawdim % 2d or 3d
case 2 % 2d
hq=quiver(real(pos),imag(pos),real(dir),imag(dir),0);
if ismember(typeofdraw,[1 3])
    ht=text(real(place),imag(place),num2str(labels));
end
case 3 % 3d
hq=quiver3(pos(:,1),pos(:,2),pos(:,3),dir(:,1),dir(:,2),dir(:,3),0);
if ismember(typeofdraw,[1 3])
    ht=text(place(:,1),place(:,2),place(:,3),num2str(labels));
end
end


if ismember(typeofdraw,[0 3])
    set(hq,'showarrowhead','off')
end

%{ 
%Example
pos =[
  -1.5543e-15 + 4.4409e-16i
   1.0000e+00 - 2.4493e-16i
   1.5000e+00 + 8.6603e-01i
   2.5000e+00 + 8.6603e-01i
   3.0000e+00 - 2.2204e-16i
   4.0000e+00 - 2.2204e-16i
   2.5000e+00 + 2.5981e+00i
   1.5000e+00 + 2.5981e+00i
   2.5000e+00 + 1.2990e+00i
   2.0000e+00 + 2.1651e+00i
   1.5000e+00 + 1.2990e+00i];
   phi=[ 2    3    4    5    6    7    8    1   11    9   10]';
 
linkdrawDez24(pos)
linkdrawDez24(pos,phi)
linkdrawDez24(pos,phi,0)
linkdrawDez24(pos,phi,1)
linkdrawDez24(pos,phi,2)
linkdrawDez24(pos,phi,3)
linkdrawDez24(pos,phi,1,0.3)
linkdrawDez24(pos,phi,1,0.3,phi)
linkdrawDez24(pos,phi,1,0.3,phi,1:8)


%}
