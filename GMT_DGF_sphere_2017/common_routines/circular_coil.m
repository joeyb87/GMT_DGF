function    [coilvec,nc,ncfull,nseg,combmat] = circular_coil(coil_radius,coil_offset,coil_rotation)
 
% Coil  definition for circular coil
% 
% Daniel Sodickson	3/23/07 
 
nc = 1; 
ncfull = 1; 
nseg = 60; 
combmat = eye(nc); 

drfac = pi/180;
phis = drfac*(0:360/nseg:(360-360/nseg));

x = coil_radius*cos(phis);
y = coil_radius*sin(phis);
z = coil_offset*ones(size(phis));

theta_coil = drfac*coil_rotation(1);
phi_coil = drfac*coil_rotation(2);
costhetacoil = cos(theta_coil);
sinthetacoil = sin(theta_coil);
cosphicoil = cos(phi_coil);
sinphicoil = sin(phi_coil);

rotation_alpha_coil = ...
    [   cosphicoil   sinphicoil   0
        -sinphicoil  cosphicoil   0
        0            0            1];
rotation_beta_coil = ...
    [   costhetacoil   0       -sinthetacoil
        0              1       0
        sinthetacoil   0       costhetacoil];

% rotation_coil = rotation_alpha_coil*rotation_beta_coil*rotation_alpha_coil';
rotation_coil = rotation_beta_coil*rotation_alpha_coil; % RL 02/28/2008 first rotate by alpha (phi) around z, then by beta (theta) around y

rotation_coil_inv = inv(rotation_coil);
r = [x;y;z];
xr = rotation_coil_inv(1,:)*r;
yr = rotation_coil_inv(2,:)*r;
zr = rotation_coil_inv(3,:)*r;
coilvec = [xr(:) yr(:) zr(:)]; 
