%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function computes the first derivative of "r" times the spherical 
%  Hankel function of "r*Ko".
%  d/dr(r*bessH(nu,r*Ko)) = bessH(nu,r*Ko) + r*d/dr(bessH(nu,r*Ko))
%                         = bessH(nu,r*Ko) + r*(Ko*(nu/r*Ko)*bessH(nu,r*Ko) - Ko*bessH(nu+1,r*Ko))
%                         = (nu + 1)*bessH(nu,r*Ko) - r*Ko*bessH(nu+1,r*Ko)
%
%  INPUT:
%  nu: order of the spherical Hankel function
%  z: argument of the spherical Hankel function
%
%  OUTPUT:
%  RderH: derivative with respect to z of the Spherical Hankel function
%
%  Name: derRspherbessH
%  Author: Riccardo Lattanzi
%  Created: 10 October 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RderH = derRspherbessH(nu, z)

RderH = -z.*spherbessH(nu + 1, z) + (nu + 1).*spherbessH(nu, z);

% END derspherbess.m
