function w = shearform(z, shear, scale, origin)
% w = shearform(z, shear, scale, origin)
% w = shearform(z, shear, rotation, translation)
% w = shearform(z, z_2, z_1, z_0)
% SHEARFORM Applies shear transformation to complex numbers, followed by
%          optional rotation and translation transformations.
%   This function takes in a set of complex numbers z and applies a shear
%   transformation to them. It also accepts rotation and translation values
%   to perform additional transformations, if desired.
%   Inputs:
%       z - A vector or matrix of complex numbers to be transformed.
%       shear          - The shear factor to be applied in the form of a complex number (it accounts for its deviation from 1i) in the sense that 1i will be mapped to shear.
%       scale       - A rotation angle in radians together with an amplification (or contraction) encoded as a complex number in polar form. If 1 or omitted, no
%                        scaling is applied.
%       translation of origin    - A complex number representing the translation to be applied.
%                        If 0 or omitted, no translation is applied.
%   Output:
%       w - The transformed complex numbers after applying shear, rotation,
%                and translation transformations.
%
% Example: if 
% z=1+1i, z0=0, z1=1, z2=1i, % then
% w=1+1i
% but if 
% z=1+1i, z0=0, z1=1, z2=1+1i, % then
% w=1i
% w=shearform(z, z_2, z_1, z_0)

z2=shear;
try  z1=scale; catch z1=1; end
try  z0=origin; catch z0=0; end

rotscale=z0+(z1-z0)*z;
sheardeviation=0; z2-(z1-z0)*1i;
m=real(sheardeviation);
n=imag(sheardeviation);
x=real(rotscale);
y=imag(rotscale);
%w=rotscale;
w=(x+m*y)+n*y*1i;