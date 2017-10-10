function plot_object_setup(s_options)

plot_1plane_flag = 1;   % 0 --> add a slider to scroll through different planes (### TO BE DONE ###)
outerfov_radius = s_options.currentradius + 0.01;

nf = 256;
np = nf;

outerfov = [2*outerfov_radius 2*outerfov_radius];
fovf = outerfov(1);
fovp = outerfov(2);

limf = -fovf/2:fovf/nf:(fovf/2 - fovf/nf);
limp = -fovp/2:fovp/np:(fovp/2 - fovp/np);

[x_fov,y_fov] = meshgrid(limf,limp);
x_fov = x_fov + ((fovf/nf)/2);
y_fov = y_fov + ((fovp/np)/2);

r_fov = sqrt(x_fov.^2 + y_fov.^2);
mask_reg4 = (r_fov <= s_options.radius3);
mask_reg3 = (r_fov <= s_options.radius2);
mask_reg2 = (r_fov <= s_options.radius1);
mask_reg1 = (r_fov > s_options.radius1);

fullfov = ones(nf,np);
% this will be used to know which voxels need to be included in the EM field calculation
fov_reg4 = fullfov.*double(mask_reg4);
fov_reg3 = fullfov.*double(mask_reg3 - mask_reg4);
fov_reg2 = fullfov.*double(mask_reg2 - mask_reg3);
fov_reg1 = fullfov.*double(mask_reg1);

multilayerfov = fov_reg1.*(fov_reg1-0.8) + fov_reg2.*(fov_reg2-0.4) + fov_reg3.*(fov_reg3-0.2) + (fov_reg4);

if plot_1plane_flag
    % Draws a circle that indicates the position of the current distribution
    % (only if its distance from the object is greater than the spatial resolution)
    if (s_options.currentradius - s_options.radius1) > s_options.radius1/nf
        region1_extension = find(multilayerfov(nf/2,(nf/2+1):end) == (fov_reg1(1)-0.8));
        current_relative_radius = round(length(region1_extension)*(s_options.currentradius - s_options.radius1)/(outerfov_radius - s_options.radius1));
        current_circle_radius = region1_extension(current_relative_radius);
        current_circle = MidpointCircle(zeros(size(multilayerfov)),current_circle_radius,nf/2,nf/2,1);
        multilayerfov = multilayerfov + current_circle;
    else
        disp('---------------------------------------------------------------------');
        disp('**WARNING** current distribution raius not shown because it would be ');
        disp('          indistinguishable from radius 1 with this voxel resolution ');
        disp('---------------------------------------------------------------------');
    end
    figure;
    set(gcf,'name','Multi-layer Sphere with Location of the Current Distribution (white circle)');
    imshow(multilayerfov,[]);
    title(['Voxel Resolution = ' sprintf('%0.2f',1000*(2*s_options.radius1)/(s_options.matrix_size(1))) ' x ' sprintf('%0.2f',1000*(2*s_options.radius1)/(s_options.matrix_size(1))) ' [mm^2]'],'FontSize',16)
    axis square
    
    figure;
    set(gcf,'name','Multi-layer FOVs');
    subplot(2,2,1)
    imshow(fov_reg4);
    title({['Region 4 (INNERMOST): r < ' sprintf('%0.2f',100*s_options.radius3) ' cm'];['Sigma = ' sprintf('%0.2f',s_options.sigma_r4) '; Eps-Rel = ' sprintf('%0.2f',s_options.epsilon_r4) ]},'FontSize',14);
    subplot(2,2,2)
    imshow(fov_reg3-0.2);
    title({['Region 3: r < ' sprintf('%0.2f',100*s_options.radius2) ' cm'];['Sigma = ' sprintf('%0.2f',s_options.sigma_r3) '; Eps-Rel = ' sprintf('%0.2f',s_options.epsilon_r3) ]},'FontSize',14);
    subplot(2,2,3)
    imshow(fov_reg2-0.4);
    title({['Region 2: r < ' sprintf('%0.2f',100*s_options.radius1) ' cm'];['Sigma = ' sprintf('%0.2f',s_options.sigma_r2) '; Eps-Rel = ' sprintf('%0.2f',s_options.epsilon_r2) ]},'FontSize',14);
    subplot(2,2,4)
    imshow(fov_reg1-0.6);
    title({['Region 1 (OUTERMOST): r > ' sprintf('%0.2f',100*s_options.radius1) ' cm'];['Sigma = 0; Eps-Rel = 1']},'FontSize',14);
    axis square
else
    disp('#### WORK IN PROGRESS ####')
end
