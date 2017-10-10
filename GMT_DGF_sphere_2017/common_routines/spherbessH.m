%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function computes the spherical Hankel function. 
%
%  INPUT:
%  nu: order of the spherical Hankel function
%  z: argument of the spherical Hankel function
%
%  OUTPUT:
%  spherbesselH: Spherical Hankel function of the first type
%
%  Name: spherbessJ
%  Author: Riccardo Lattanzi
%  Created: 10 October 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spherbesselH = spherbessH(nu, z)

spherbesselH = sqrt(pi/2) .* sqrt(1./z) .* besselh((nu + 1/2), 1, z);
