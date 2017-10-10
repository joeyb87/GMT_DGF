%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function computes the first derivative of "r" times the spherical 
%  Bessel function of "r*Ko".
%  d/dr(r*bessJ(nu,r*Ko)) = bessJ(nu,r*Ko) + r*d/dr(bessJ(nu,r*Ko))
%                         = bessJ(nu,r*Ko) + r*(Ko*(nu/r*Ko)*bessJ(nu,r*Ko) - Ko*bessJ(nu+1,r*Ko))
%                         = (nu + 1)*bessJ(nu,r*Ko) - r*Ko*bessJ(nu+1,r*Ko)
%
%  INPUT:
%  nu: order of the spherical Bessel function
%  z: argument of the spherical Bessel function
%
%  OUTPUT:
%  RderJ: derivative with respect to z of the Spherical Bessel function
%
%  Name: derRspherbessJ
%  Author: Riccardo Lattanzi
%  Created: 27 March 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RderJ = derRspherbessJ(nu, z)

RderJ = -z.*spherbessJ(nu + 1, z) + (nu + 1).*spherbessJ(nu, z);

% END derspherbess.m
