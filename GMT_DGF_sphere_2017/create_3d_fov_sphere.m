function [x_fov,y_fov,z_fov,mask_reg1,mask_reg2,mask_reg3,mask_reg4] = create_3d_fov_sphere(s_options)

fovx = 2*s_options.radius1;
fovy = fovx;
fovz = fovx;

nx = s_options.matrix_size(1);
ny = s_options.matrix_size(2);
nz = s_options.matrix_size(3);

gridlimx = -fovx/2:fovx/nx:(fovx/2 - fovx/nx);
gridlimy = -fovy/2:fovy/ny:(fovy/2 - fovy/ny);
gridlimz = -fovz/2:fovz/nz:(fovz/2 - fovz/nz);

[x_fov,y_fov,z_fov] = meshgrid(gridlimx,gridlimy,gridlimz);
x_fov = x_fov + ((fovx/nx)/2);
y_fov = y_fov + ((fovy/ny)/2);
z_fov = z_fov + ((fovz/nz)/2);

r_fov = sqrt(x_fov.^2 + y_fov.^2 + z_fov.^2);
mask_reg4 = (r_fov <= s_options.radius3);
mask_reg3 = (r_fov <= s_options.radius2) - (r_fov <= s_options.radius3);
mask_reg3 = logical(mask_reg3);
mask_reg2 = (r_fov <= s_options.radius1) - (r_fov <= s_options.radius2);
mask_reg2 = logical(mask_reg2);
mask_reg1 = (r_fov > s_options.radius1);

