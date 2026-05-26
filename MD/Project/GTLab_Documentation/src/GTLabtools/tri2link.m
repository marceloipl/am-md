function [g, phi, s] = tri2link(T, V)
% [g, phi, s] = tri2link(T, V)
%
% Wed, 03Oct2024 --- small changes relative to the previous version of 15Jun2024
%
%TRI2LINK Converts a mesh (T,V) into an anchored link (g,phi,s).
%   [g, phi, s] = tri2link(T, V) takes a mesh represented by triangle
%   indices T and vertex coordinates V and converts it into an anchored link
%   structure with geometrical information g, edge-connecting map phi, and
%   starting vertex map s.
%
%   Inputs:
%       T: Matrix representing triangular faces of the mesh. Each row
%          corresponds to a triangle, and each column contains indices
%          referencing vertex coordinates in V.
%       V: Matrix containing vertex coordinates of the mesh, where each row
%          represents a vertex.
%
%   Outputs:
%       g: Vector containing geometrical information for each edge.
%          Calculated as g(a) = log(r(a)) + i*t(a), where r(a) is the ratio
%          of edge lengths and t(a) is the difference in edge angles.
%       phi: Vector representing the map connecting edges along the same face.
%       s: Vector representing the map assigning starting vertices to each edge.
%
%   The function processes the input mesh and constructs the anchored link
%   structure using the provided vertex coordinates and triangle definitions.
%   This structure allows for further analysis and manipulation of the mesh in
%   the context of the anchored link framework.
% Wed, 03Oct2024, martins.ferreira@ipleiria.pt


id=@(x) (1:numel(x))';    % identity map with the same domain as vector x
[nT,nf]=size(T);   % nT=number of triangles (or faces); nf= number of edges in each face (nf=3 for triangles)

s=T'; % T=[a(:) b(:) c(:)]
%t=s([2:nf 1],:); % t can be recovered as t=s(phi);
s=s(:); %t=t(:); % working with source and target as column vectors
%idA=id(s); % the identity on edges
phi=reshape(1:nf*nT,nf,nT);
phi=phi([2:nf 1],:);
phi=phi(:); % phi is such that t=s(phi)
%source=s;  % so not to overuse the letter s

% computing g, based on the amplitwist approach but making the necessary changes
%source=s;
t=s(phi); idA=id(s);
dir=@(a) V(t(a),:)-V(s(a),:); % direction of vector associated with edge indexed by a, interpreted as a:s(a)-->t(a)
norm=@(v) sqrt(dot(v',v')');   % v must be of size n-by-3
ampli=norm(dir(phi))./norm(dir(idA));
twist=acos(dot(dir(phi)',dir(idA)')'./(norm(dir(phi)).*norm(dir(idA))));
%g=ampli.*exp(1i*twist);  %older interpretation
g=log(ampli)+1i*twist;

end % function tri2link

%{
Example
[V,T] = platonic_solids(2);  % dodecahedron --- available from Matlab Exchange files
[g,phi,s]=tri2link(T,V)

%}
