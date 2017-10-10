% distribute_coils.m, 30.6,4 | Florian Wiesinger
%
% function [theta,phi]=distribute_coils(N)
%
% This function distributes n identically charged particles evenly over 
% the unit-sphere. This is achieved by the minimization of the 
% corresponding e.g. electrostatic potential sum_inv_rij. This is used
% for the uniform distribution of N coils.
%
% INPUT:
% N         ... number of particles / coils
%
% OUTPUT:
% phi       ... PHI   - coordinates of evenly distributed particles / coils
% theta     ... THETA - coordinates of evenly distributed particles / coils

function [theta,phi]=distribute_coils(N)

% MINIMIZATION OF e.g. an electrostatic POTENTIAL of particles (/coils) placed on a unit sphere: sum{(1/rij)^2}
x0=pi*[2*rand(N-1,1),rand(N-1,1)];      % phi and theta angles for initializing the minimization
[x,resnorm,residual,exitflag] = lsqnonlin(@sum_inv_rij,x0);
if exitflag <= 0
    display('!! WARNING: DISTRIBUTE COILS COULD NOT FIND A MINIMUM !!');
end
phi=[0;x(:,1)]; theta=[pi/2;pi/2;x(2:end,2)]; 

function y=sum_inv_rij(x)
% FIRST PARTICLE / COIL IS FIXED AT: phi=0, theta=pi/2 !!
phi=[0;squeeze(x(:,1))]; deta=[pi/2;pi/2;squeeze(x(2:end,2))]; n=length(phi);
x=sin(deta).*cos(phi); y=sin(deta).*sin(phi); z=cos(deta);
[xm1,xm2]=meshgrid(x,x); [ym1,ym2]=meshgrid(y,y); [zm1,zm2]=meshgrid(z,z);
y = reshape((ones(n) - eye(n)) ./ ((xm1-xm2).^2 + (ym1-ym2).^2 + (zm1-zm2).^2 + eye(n)),1,n^2);

