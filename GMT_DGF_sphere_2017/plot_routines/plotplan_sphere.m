function	plotplan_sphere(x,y,z)

% Plots plane of interest as a patch on current axis defined by the edges
% of coordinate arrays x,y,z.
% Used by coilsim.m.
%
% function	plotplan(x,y,z)
% 
% input:
%	x,y,z: coordinate arrays defining plane of interest
%
% Daniel Sodickson
% Version History:
% 1.2: 10/13/98
% 1.1: 8/8/98

[nz,ny] = size(x);

hp = patch(	[x(1,1); x(1,ny); x(nz,ny); x(nz,1)],...
		[y(1,1); y(1,ny); y(nz,ny); y(nz,1)],...
		[z(1,1); z(1,ny); z(nz,ny); z(nz,1)],...
		[z(1,1); z(1,ny); z(nz,ny); z(nz,1)]	);
% set(hp,'erasemode','normal','facecolor','flat')

set(hp,'facecolor',[0.4 0.4 0.4])
% set(hp,'erasemode','normal','facecolor',[0.4 0.4 0.4])

% set(hp,'erasemode','normal','facecolor',[1 1 0]) % for plots used in ISMRM talk
% set(hp,'erasemode','normal','facecolor',[0.8 0.94 0.98]) % for plots used in ISMRM talk


