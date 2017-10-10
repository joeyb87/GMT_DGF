        
function  [output_values] = dgf_sphere_calc_snr(... 
                  whichcurrents,fieldstrength,acceleration,expsnr_num,s_opts,g_opts,path_opts,currentpatternmatrixsize,fieldpatternmatrixsize,...
                  fieldpatternthreeplaneoffset,whichvoxels,whichvoxelweights,plot_poynting_sphere_flag,plot_regions_fov)

%--------------------------------------------------------------------------
% 
% function  [output_values] = dgf_sphere_calc_snr(... 
%                  whichcurrents,fieldstrength,acceleration,simulation_options,geometry_options,path_options,currentpatternmatrixsize,whichvoxels,whichvoxelweights);
% ------
% input:
% ------
%   whichcurrents: choice of current basis functions to include (1 = divergence-free mag dipole, 2 = curl-free electric dipole) [default [1 2] = both]
%   fieldstrength: field strength (Tesla)
%   acceleration: 1 x 2 vector containing acceleration factors in each of two dimensions
%   expsnr_num: experimental SNR scaling
%   s_opts:
%       lmax: 2*[(lmax + 1)^2 - 1] modes in the harmonic expansion (NB: l=0 is not included because it gives no SNR contribution)
%       include_ult_conductor_losses: 1 --> include losses due to Johnson noise (by uSNR definition should be 0), 0 --> sigma infinite and no coil noise
%       include_coil_conductor_losses: 1 --> include coil losses due to Johnson noise, 0 --> sigma infinite and no coil noise
%       include_ult_radiation_losses: 1 --> include modes' radiation losses (average Poynting vector, no shielding)
%       include_coil_radiation_losses: 1 --> include coils' radiation losses (average Poynting vector, no shielding)
%       include_circuit_losses: 1 --> include noise (not for ultimate case) due to impedance transformation network and LNA (** currently only noise figure **)
%       include_psi_coil: 1 --> include noise correlation matrix in coil optimization
%       include_experimental_snr_scaling: 1 --> include sequence related scaling factors
%       noisefactor: noise factor of the receive chain
%       compute_ultSNR_flag: 1 --> compute ultimate intrinsic SNR
%       compute_coilSNR_flag: 1 --> calculate SNR and g also for circular surface coils
%       compute_ult_current_pattern_flag: 1 --> compute ultimate current pattern
%       compute_coil_current_pattern_flag: 1 --> compute coil current pattern
%       save_efields_ult_flag: 1 --> save a matrix with the modes' electric fields at each position
%       save_efields_coil_flag: 1 --> save a matrix with the coils' electric fields at each position
%       save_bfields_ult_flag: 1 --> save a matrix with the modes' magnetic fields at each position
%       save_bfields_coil_flag: 1 --> save a matrix with the coils' magnetic fields at each position
%       save_ult_current_weights_flag: 1 --> save ultimate current weights
%       save_coil_current_weights_flag: 1 --> save circular surface coil current weights
%       save_ult_currents_flag: 1 --> save ultimate current patterns
%       save_coil_currents_flag: 1 --> save coil current patterns
%       save_noise_contributions: 1 --> save sample, coil and radiation losses separately 
%       sigma_coil: conductivity of coil material, for estimates of coil noise (Siemens/m)
%       d_coil: thickness of coil material, for estimates of coil noise (m)
%       tissue_type: 1 --> brain from Wiesinger et al, 2 --> dog skeletal muscle Cole-Cole-scaled from Stoy et al [default 1], 3 --> Dog skeletal muscle from Schnell thesis Appendix C.
%       epsilon_rel: relative permittivity of imaged body  [NaN --> select tissue values based on field strength]
%       sigma: conductivity of imaged body (Siemens/m) [NaN --> select tissue values based on field strength]


%   g_opts:
%       preset_fov_geometry_flag: 1 --> load FOV data from .MAT file
%       fov_file: file with FOV data (if [] user is allow to browse directories)
%	    fovf: fov in the frequency-encoding direction
%	    fovp: fov in the phase-encoding direction
%       matrix_size: 1 x 2 vector containing image matrix dimensions
%       radius3: radius of dielectric sphere (m)
%       currentradius: radius of coil former (m)
%	    phasedir: phase encoding / foldover direction flag ('FH' or 'LR')
%	    patientposition: 'headfirst' or 'feetfirst'
%	    patientorientation: 'supine','prone','ldecub', or 'rdecub'
%	    sliceorientation: reference slice orientation ('transverse','sagittal', or 'coronal')
%       image_plane_offset: [apoff lroff fhoff] list of linear offsets in the AP, LR, FH direction
%       image_plane_orientation: [apang lrang fhang] list of angular offsets about the AP, LR, FH axis
%       ncoils: number of circular surface coils in the array
%       coil_rotations: ncoils x 2 array of inclination and azimuth [theta=beta phi=alpha] of coil center (degrees)
%       coil_radii: ncoils x 1 array of radii of circular coils (m)
%       coil_offsets: ncoils x 1 array of displacements of coil center from sphere center (m)
%   path_opts:
%       commondir: directory with general routines used in the calculations
%       plotdir: directory with plotting routines
%       basissetdir: sub-directory to store results for the ultimate case
%       circcoildir: sub-directory to store results for the coil case
%       moviedir: sub-directory to store movies
%       logdir: sub-directory to store log files
%       geometrydir: sub-directory with preset coil geometries
%   currentpatternmatrixsize: dimensions of stored current pattern for each voxel position [default 32 x 32]                
%   whichvoxels: (number of chosen voxels x 2) array of voxel indices for which to compute current patterns [default NaN --> save current patterns for all voxels]
%   whichvoxelweights: (number of chosen voxels x 2) array of voxel indices for which to save weights [default NaN --> save weights for all voxels]
%   plot_poynting_sphere
%
% -------
% output:
% -------
%   output_values:
%       snr_ult: matrix_size(1) x matrix_size(2) array of SNR values
%       g_ult: matrix_size(1) x matrix_size(2) array of g-factor values
%       snr_coil: matrix_size(1) x matrix_size(2) array of circular surface coil array SNR values
%       g_coil: matrix_size(1) x matrix_size(2) array of circular surface coil array g-factor values
%       psi_coil: coil noise covariance matrix
%       mask: matrix_size(1) x matrix_size(2) binary mask showing object boundary
%       whichvoxels: (number of chosen voxels x 2) array of voxel indices for which to compute current patterns [default NaN --> save current patterns for all voxels]
%       whichvoxelweights: (number of chosen voxels x 2) array of voxel indices for which to save weights [default NaN --> save weights for all voxels]
%       weights: length(whichn) x length(whichm) x length(whichcurrents) x size(whichvoxelweights) array of optimal weights for current basis elements
%       weights_coil: length(whichn) x length(whichm) x length(whichcurrents) x size(whichvoxelweights) array of optimal weights for circular surface coils
%       currentpattern: currentpatternmatrixsize(1) x currentpatternmatrixsize(2) x 2 x number of chosen voxels array of current pattern values as a function of phi and theta
%       currentpattern_coil: currentpatternmatrixsize(1) x currentpatternmatrixsize(2) x 2 x number of chosen voxels array of circular surface current pattern values as a function of phi and z
%       currentphi: phi coordinates for current pattern plotting
%       currenttheta: theta coordinates for current pattern plotting
%       epsilon_rel: the value of epsilon_rel actually used if epsilon_rel was set to NaN on input
%       sigma: the value of sigma actually used if sigma was set to NaN on input
%       x_fov: x coordinates of chosen fov
%       y_fov: y coordinates of chosen fov
%       z_fov: z coordinates of chosen fov
%
% Riccardo Lattanzi, June 1, 2012
%
% Debugging notes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(path_opts.commondir);
addpath(path_opts.plotdir);
addpath(path_opts.tissuedir);
addpath(path_opts.fovdir);

matrix_size = g_opts.matrix_size;

lmax = s_opts.lmax;
snr_radius = s_opts.snr_radius;
currentradius = g_opts.currentradius;
radius3 = g_opts.radius3;
radius2 = g_opts.radius2;
radius1 = g_opts.radius1;
sigma_coil = s_opts.sigma_coil;
ncoils = g_opts.ncoils;
coil_rotations = g_opts.coil_rotations;
coil_radii = g_opts.coil_radii;
coil_offsets = g_opts.coil_offsets;
tissue_r2 = s_opts.tissue_r2;
tissue_r3 = s_opts.tissue_r3;
tissue_r4 = s_opts.tissue_r4;

% RL 9/26/2013 to account for the 90 degrees shift between UISNR
% and coil SNR
if s_opts.compute_coilSNR_flag
    shiftrotation = zeros(size(coil_rotations));
    shiftrotation(:,2) = 90*ones(size(coil_rotations,1),1);
    coil_rotations = coil_rotations - shiftrotation;
end
%%--------------------%%

disp('***   running dgf_sphere_calc_snr...')
tic

if all(whichcurrents==1),
    disp('***   Magnetic dipole currents only');
elseif all(whichcurrents==2),
    disp('***   Electric dipole currents only');
else
    disp('***   Magnetic and electric dipole currents');    
end
disp('***---------------------------------------------***');

if ~s_opts.compute_coilSNR_flag
    PRdivPL = [];
end

%---------------------------------%
%--- DEFINE FOV OF ALL REGIONS ---%
%---------------------------------%

if g_opts.preset_fov_geometry_flag
    load(g_opts.fov_file);
    nf = matrix_size(1);
    np = matrix_size(2);
    maxrad = sqrt(max(x_fov(:)).^2 + max(y_fov(:)).^2 + max(z_fov(:)).^2);
    fovf = max(x_fov(:)); % OK only for transverse: need to be fixed!
    fovp = max(y_fov(:));
    if maxrad < g_opts.radius3
        s_opts.save_emfields_ult_r3_flag = 0;
        s_opts.save_emfields_coil_r3_flag = 0;
        s_opts.save_emfields_ult_r2_flag = 0;
        s_opts.save_emfields_coil_r2_flag = 0;
        s_opts.save_emfields_ult_r1_flag = 0;
        s_opts.save_emfields_coil_r1_flag = 0;
    end
    if maxrad < g_opts.radius2
        s_opts.save_emfields_ult_r2_flag = 0;
        s_opts.save_emfields_coil_r2_flag = 0;
        s_opts.save_emfields_ult_r1_flag = 0;
        s_opts.save_emfields_coil_r1_flag = 0;
    end
    if maxrad < g_opts.radius1
        s_opts.save_emfields_ult_r1_flag = 0;
        s_opts.save_emfields_coil_r1_flag = 0;
    end
else
    outerfov = [g_opts.mask_radius*2*g_opts.outerfov_radius g_opts.mask_radius*2*g_opts.outerfov_radius];
    fovf = outerfov(1);
    fovp = outerfov(2);
    nf = matrix_size(1);
    np = matrix_size(2);
    apoff = g_opts.image_plane_offset(1); lroff = g_opts.image_plane_offset(2); fhoff = g_opts.image_plane_offset(3);
    apang = g_opts.image_plane_orientation(1); lrang = g_opts.image_plane_orientation(2); fhang = g_opts.image_plane_orientation(3);
    
    [x_fov,y_fov,z_fov] = sphere_mkplane(fovf,fovp,nf,np,g_opts.phasedir,...
        g_opts.patientposition,g_opts.patientorientation,...
        g_opts.sliceorientation,apoff,lroff,fhoff,apang,lrang,fhang);
    
    % Make pixels positions symmetric with respect to the origin for a transverse plane
    switch g_opts.sliceorientation
        case 'coronal'
            if g_opts.phasedir == 'FH'
                x_fov = x_fov + ((fovf/nf)/2);
                z_fov = z_fov + ((fovp/np)/2);
            elseif g_opts.phasedir == 'LR'
                x_fov = x_fov + ((fovp/np)/2);
                z_fov = z_fov + ((fovf/nf)/2);
            end
        case 'transverse'
            if g_opts.phasedir == 'FH'
                x_fov = x_fov + ((fovf/nf)/2);
                y_fov = y_fov + ((fovp/np)/2);
            elseif g_opts.phasedir == 'LR'
                x_fov = x_fov + ((fovp/np)/2);
                y_fov = y_fov + ((fovf/nf)/2);
            end
        case 'sagittal'
            if g_opts.phasedir == 'FH'
                y_fov = y_fov + ((fovf/nf)/2);
                z_fov = z_fov + ((fovp/np)/2);
            elseif g_opts.phasedir == 'LR'
                y_fov = y_fov + ((fovp/np)/2);
                z_fov = z_fov + ((fovf/nf)/2);
            end
    end
end

%     % TEST
%     y_fov = z_fov + 0.0015;

% % Define a mask for the pixel within the circular section
% r_fov = sqrt(x_fov.^2 + y_fov.^2);
% mask = (r_fov <= radius3);

% Define a mask for the pixel within the circular section
switch  g_opts.sliceorientation
    case 'coronal'
        r_fov = sqrt(x_fov.^2 + z_fov.^2);
    case 'transverse'
        r_fov = sqrt(x_fov.^2 + y_fov.^2);
    case 'sagittal'
        r_fov = sqrt(y_fov.^2 + z_fov.^2);
end
mask_reg4 = (r_fov <= radius3);
mask_reg3 = (r_fov <= radius2);
mask_reg2 = (r_fov <= radius1);
mask_reg1 = (r_fov > radius1);

fullfov = ones(matrix_size);
% this will be used to know which voxels need to be included in the EM field calculation
fov_reg4 = fullfov.*double(mask_reg4);
fov_reg3 = fullfov.*double(mask_reg3 - mask_reg4);
fov_reg2 = fullfov.*double(mask_reg2 - mask_reg3);
fov_reg1 = fullfov.*double(mask_reg1);
switch snr_radius
    % this will be used to know which voxels need to be included in the SNR calculation
    case 1
        snr_fov = fullfov.*double(mask_reg2);
        mask_snr = mask_reg2;
        disp('')
        disp('-- SNR will be calculate for Region 4, 3 and 2 --')
    case 2
        snr_fov = fullfov.*double(mask_reg3);
        mask_snr = mask_reg3;
        disp('')
        disp('-- SNR will be calculate for Region 4 and 3 --')
    case 3
        snr_fov = fov_reg4;
        mask_snr = mask_reg4;
        disp('')
        disp('-- SNR will be calculate only for Region 4 --')
end

%-- check if all varaiables are there --%
if isempty(x_fov) || isempty(y_fov) || isempty(z_fov)
    disp('-------------------------------------------------------');
    disp('**ERROR** unspecified FOV variable');
    disp('-------------------------------------------------------');
    pause
end

%------------------------------------------%
%--- TISSUE PROPERTIES AND WAVE VECTORS ---%
%------------------------------------------%

mu = 4*pi*1e-7;         % permeability of free space [Wb][A^-1][m^-1]
c = 3e8;                % speed of light [m][s]
epsilon_0 = 1/(mu*c^2); % permittivity [C][V^-1] = [s^4][A^2][m^-2][kg^-2]

omega = 2*pi*42.576e6*fieldstrength; % Larmor frequency [Hz]

if isnan(s_opts.sigma_r4) || isnan(s_opts.epsilon_r4),
    tissuefile_r4 = [tissue_r4 '.mat'];
    load(tissuefile_r4);
    fieldset = tissueproperties(1,:)./42.576e6;
    sigma_set = tissueproperties(2,:);
    epsilon_rel_set = tissueproperties(3,:);
    epsilon_r4 = spline(fieldset,epsilon_rel_set,fieldstrength);
    sigma_r4 = spline(fieldset,sigma_set,fieldstrength);
    clear tissueproperties
else
    epsilon_r4 = s_opts.epsilon_r4;
    sigma_r4 = s_opts.sigma_r4;
end
if isnan(s_opts.sigma_r3) || isnan(s_opts.epsilon_r3),
    tissuefile_r3 = [tissue_r3 '.mat'];
    load(tissuefile_r3);
    fieldset = tissueproperties(1,:)./42.576e6;
    sigma_set = tissueproperties(2,:);
    epsilon_rel_set = tissueproperties(3,:);
    epsilon_r3 = spline(fieldset,epsilon_rel_set,fieldstrength);
    sigma_r3 = spline(fieldset,sigma_set,fieldstrength);
    clear tissueproperties
else
    epsilon_r3 = s_opts.epsilon_r3;
    sigma_r3 = s_opts.sigma_r3;
end
if isnan(s_opts.sigma_r2) || isnan(s_opts.epsilon_r2),
    tissuefile_r2 = [tissue_r2 '.mat'];
    load(tissuefile_r2);
    fieldset = tissueproperties(1,:)./42.576e6;
    sigma_set = tissueproperties(2,:);
    epsilon_rel_set = tissueproperties(3,:);
    epsilon_r2 = spline(fieldset,epsilon_rel_set,fieldstrength);
    sigma_r2 = spline(fieldset,sigma_set,fieldstrength);
    clear tissueproperties
else
    epsilon_r2 = s_opts.epsilon_r2;
    sigma_r2 = s_opts.sigma_r2;
end

% -- create results file name -- 
if length(whichcurrents) > 1
    curr_label = 3;
else
    curr_label = whichcurrents;
end
filetosave_ult = [ path_opts.basissetdir 'ultimate_snr_results/UISNR_sphere_r[' ...
    num2str(radius3*100) '-' num2str(g_opts.radius2*100) '-' num2str(g_opts.radius1*100) '-' num2str(currentradius*100) 'cm]_' ...
    sprintf('%2.1fT_',fieldstrength) num2str(matrix_size(1)) 'x' num2str(matrix_size(2)) ...
    '_acc_' num2str(acceleration(1)) 'x' num2str(acceleration(2)) ...
    '_e[' sprintf('%1.0f-',epsilon_r4) sprintf('%1.0f-',epsilon_r3) sprintf('%1.0f',epsilon_r2) ']' ...
    '_s[' sprintf('%2.2f-',sigma_r4) sprintf('%2.2f-',sigma_r3) sprintf('%2.2f',sigma_r2) ']' ...
    '_l' num2str(lmax) '_cur_' num2str(curr_label) '.mat'];
% for loop coils curr_label is always equal to one, so it's not included in the file name
if s_opts.compute_coilSNR_flag
    if s_opts.include_psi_coil
        filetosave_coil = [ path_opts.circcoildir 'coils_snr_results/cSNR_sphere_r[' ...
            num2str(radius3*100) '-' num2str(g_opts.radius2*100) '-' num2str(g_opts.radius1*100) '-' num2str(currentradius*100) 'cm]_' ...
            sprintf('%2.1fT_',fieldstrength) num2str(matrix_size(1)) 'x' num2str(matrix_size(2)) ...
            '_acc_' num2str(acceleration(1)) 'x' num2str(acceleration(2)) ...
            '_e[' sprintf('%1.0f-',epsilon_r4) sprintf('%1.0f-',epsilon_r3) sprintf('%1.0f',epsilon_r2) ']' ...
            '_s[' sprintf('%2.2f-',sigma_r4) sprintf('%2.2f-',sigma_r3) sprintf('%2.2f',sigma_r2) ']' ...
            '_c' num2str(ncoils) '_l' num2str(lmax) '.mat'];
    else
        filetosave_coil = [ path_opts.circcoildir 'coils_snr_results/cSNR_sphere_r['...
            num2str(radius3*100) '-' num2str(g_opts.radius2*100) '-' num2str(g_opts.radius1*100) '-' num2str(currentradius*100) 'cm]_' ...
            sprintf('%2.1fT_',fieldstrength) num2str(matrix_size(1)) 'x' num2str(matrix_size(2)) ...
            '_acc_' num2str(acceleration(1)) 'x' num2str(acceleration(2)) ...
            '_e[' sprintf('%1.0f-',epsilon_r4) sprintf('%1.0f-',epsilon_r3) sprintf('%1.0f',epsilon_r2) ']' ...
            '_s[' sprintf('%2.2f-',sigma_r4) sprintf('%2.2f-',sigma_r3) sprintf('%2.2f',sigma_r2) ']' ...
            '_c' num2str(ncoils) '_l' num2str(lmax) '_nopsi.mat'];
    end
end

epsilon_r4 = epsilon_r4*epsilon_0;
epsilon_r3 = epsilon_r3*epsilon_0;
epsilon_r2 = epsilon_r2*epsilon_0;

disp('***---------------------------------------------***');
disp(['B_o = ' num2str(fieldstrength) ' [T]']);
disp(['omega = ' num2str(omega/1E6) ' [MHz]']);
disp(['acceleration factor = ' num2str(acceleration(1)) 'x' num2str(acceleration(2)) ]);
if size(x_fov,1) > 1
    disp(['Matrix Size = ' num2str(nf) ' x '  num2str(np)]);
    disp(['Voxel Resolution = ' sprintf('%0.2f',1000*fovf/nf) ' x ' sprintf('%0.2f',1000*fovp/np) ' [mm^2]']);
end
disp('REGION 4 (Innermost)');
disp(['       radius 4: 0 < r < ' num2str(radius3*100,3) ' [cm]']);
disp(['       sigma 4 = ' num2str(sigma_r4) ' [ohm^-1][m^-1]']);
disp(['       epsilon rel 4 = ' num2str(epsilon_r4/epsilon_0)]);
disp('REGION 3');
disp(['       radius 3: ' num2str(radius3*100,3) ' < r < ' num2str(radius2*100,3) ' [cm]']);
disp(['       sigma 3 = ' num2str(sigma_r3) ' [ohm^-1][m^-1]']);
disp(['       epsilon rel 3 = ' num2str(epsilon_r3/epsilon_0)]);
disp('REGION 2');
disp(['       radius 2: ' num2str(radius2*100,3) ' < r < ' num2str(radius1*100,3) ' [cm]']);
disp(['       sigma 2 = ' num2str(sigma_r2) ' [ohm^-1][m^-1]']);
disp(['       epsilon rel 2 = ' num2str(epsilon_r2/epsilon_0)]);
disp('REGION 1 (Outermost)');
disp(['       radius 1: r > ' num2str(radius1*100,3) ' [m]']);
disp([' with current distribution at r = ' num2str(currentradius*100,3) ' [cm]']);
disp(['       sigma 1 = 0']);
disp(['       epsilon rel 1 = 1']);
disp('***---------------------------------------------***');

if plot_regions_fov
    multilayerfov = fov_reg1.*(fov_reg1-0.8) + fov_reg2.*(fov_reg2-0.4) + fov_reg3.*(fov_reg3-0.2) + (fov_reg4);
    % Draws a circle that indicates the position of the current distribution
    % (only if its distance from the object is greater than the spatial resolution)
    if (currentradius - radius1) > fovf/nf 
        region1_extension = find(multilayerfov(nf/2,(nf/2+1):end) == (fov_reg1(1)-0.8));
        current_relative_radius = round(length(region1_extension)*(currentradius - radius1)/(g_opts.outerfov_radius - radius1));
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
    title(['Voxel Resolution = ' sprintf('%0.2f',1000*fovf/nf) ' x ' sprintf('%0.2f',1000*fovp/np) ' [mm^2]'],'FontSize',16)
    axis square
    
    figure;
    set(gcf,'name','Multi-layer FOVs');
    subplot(2,2,1)
    imshow(fov_reg4);
    title({['Region 4 (INNERMOST): r < ' sprintf('%0.2f',100*radius3) ' cm'];['Sigma = ' sprintf('%0.2f',sigma_r4) '; Eps-Rel = ' sprintf('%0.2f',epsilon_r4/epsilon_0) ]},'FontSize',14);
    subplot(2,2,2)
    imshow(fov_reg3-0.2);
    title({['Region 3: r < ' sprintf('%0.2f',100*radius2) ' cm'];['Sigma = ' sprintf('%0.2f',sigma_r3) '; Eps-Rel = ' sprintf('%0.2f',epsilon_r3/epsilon_0) ]},'FontSize',14);
    subplot(2,2,3)
    imshow(fov_reg2-0.4);
    title({['Region 2: r < ' sprintf('%0.2f',100*radius1) ' cm'];['Sigma = ' sprintf('%0.2f',sigma_r2) '; Eps-Rel = ' sprintf('%0.2f',epsilon_r2/epsilon_0) ]},'FontSize',14);
    subplot(2,2,4)
    imshow(fov_reg1-0.6);
    title({['Region 1 (OUTERMOST): r > ' sprintf('%0.2f',100*radius1) ' cm'];['Sigma = 0; Eps-Rel = 1']},'FontSize',14);
    axis square
end

k_0_squared = omega*omega*epsilon_0*mu;
% k_0_squared = omega*omega*500*epsilon_0*mu;
k_0 = sqrt(k_0_squared);

k_4_squared = omega*mu*(omega*epsilon_r4+1i*sigma_r4);
k_4 = sqrt(k_4_squared);

k_3_squared = omega*mu*(omega*epsilon_r3+1i*sigma_r3);
% k_3_squared = omega*omega*epsilon_0*mu;
k_3 = sqrt(k_3_squared);

k_2_squared = omega*mu*(omega*epsilon_r2+1i*sigma_r2);
% k_2_squared = omega*omega*epsilon_0*mu;
k_2 = sqrt(k_2_squared);

%--------------------------------------------
%        Variables Initialization
%--------------------------------------------

if all(whichcurrents == 1)  || all(whichcurrents == 2),
    numbasis = (lmax + 1)^2 - 1;
    counterincrement = 0;
    phi_matrix_set = zeros(1,numbasis);
    P_L_matrix_set  = zeros(2,2*numbasis); % P_L is always 2x2
    if s_opts.include_ult_conductor_losses || s_opts.include_coil_conductor_losses
        P_A_matrix_set  = zeros(numbasis,1);
    end
    if s_opts.include_ult_radiation_losses || s_opts.include_coil_radiation_losses
        P_R_matrix_set  = zeros(numbasis,1);
    end
else
    numbasis = 2*((lmax + 1)^2 - 1);
    counterincrement = 1;
    phi_matrix_set = zeros(2,numbasis);
    P_L_matrix_set  = zeros(2,numbasis);
    if s_opts.include_ult_conductor_losses || s_opts.include_coil_conductor_losses
        P_A_matrix_set  = zeros(numbasis,1);
    end
    if s_opts.include_ult_radiation_losses || s_opts.include_coil_radiation_losses
        P_R_matrix_set  = zeros(numbasis,1);
    end
%     if compute_ult_magn_and_ele
%         halfnumbasis = nn*nm;
%     end
end

if s_opts.compute_coilSNR_flag
    drfac = pi/180;
    rot_coil = coil_rotations*drfac;
end

if isnan(s_opts.d_coil),
    d_coil = 1/sqrt(omega*mu*sigma_coil/2);
    % use skin depth of coil conductor as coil thickness
    % d_coil = delta_coil = 1/sqrt(pi*f*mu_0*sigma_coil);
else
    d_coil = s_opts.d_coil;
end
if s_opts.save_noise_contributions
    P_L_set_r4 = zeros(numbasis,1);
    P_L_set_r3 = zeros(numbasis,1);
    P_L_set_r2 = zeros(numbasis,1);
    P_A_set = zeros(numbasis,1);
    P_R_set = zeros(numbasis,1);
else
    P_L_set_r4 = [];
    P_L_set_r3 = [];
    P_L_set_r2 = [];
    P_A_set = [];
    P_R_set = [];
end

% -- set up for saving current patterns --
if s_opts.compute_ult_current_pattern_flag
    if isempty(whichvoxels),
        [ind1,ind2] = meshgrid(1:matrix_size(1),1:matrix_size(2));
        whichvoxels = [ind1(:) ind2(:)];
    end
    currentpattern = zeros([currentpatternmatrixsize 2 size(whichvoxels,1)]);
    if s_opts.compute_coil_current_pattern_flag
        currentpattern_coil = currentpattern;
    else
        currentpattern_coil = [];
    end
    delta_phi = 2*pi/currentpatternmatrixsize(1);
    delta_theta = pi/currentpatternmatrixsize(2);                                                 
%     [currentphi,currenttheta] = ndgrid(-pi:delta_phi:(pi-delta_phi),-pi/2:delta_theta:(pi/2-delta_theta));
    [currentphi,currenttheta] = ndgrid(0:delta_phi:(2*pi-delta_phi),0:delta_theta:(pi-delta_theta));
    cos_currenttheta = cos(currenttheta);
%     cot_currenttheta = cot(currenttheta);
    csc_currenttheta = csc(currenttheta); % NB: it's infinite at theta = 0
else
    currentpattern = [];
    currentpattern_coil = [];
    currentphi = [];
    currenttheta = [];
end
% -- set up for saving optimal voxel weights --
if s_opts.save_ult_current_weights_flag == 1,
    if isempty(whichvoxelweights),
        [ind1,ind2] = meshgrid(1:matrix_size(1),1:matrix_size(2));
        whichvoxelweights = [ind1(:) ind2(:)];
    end
    weights = zeros([((lmax+1)^2 - 1) length(whichcurrents) size(whichvoxelweights,1)]);
else
    weights = [];
end
if s_opts.save_coil_current_weights_flag == 1,
    weights_coil = zeros([ncoils size(whichvoxelweights,1)]);
else
    weights_coil = [];
end

% -- set up for saving field patterns --
if s_opts.compute_ult_field_patterns || s_opts.compute_coil_field_patterns
    if isempty(whichvoxels), % --> save field patterns for all voxels
        [ind1,ind2] = meshgrid(1:matrix_size(1),1:matrix_size(2));
        whichvoxels = [ind1(:) ind2(:)];
    end
    fpms = fieldpatternmatrixsize;
    fieldxoffset = fieldpatternthreeplaneoffset(1);
    fieldyoffset = fieldpatternthreeplaneoffset(2);
    fieldzoffset = fieldpatternthreeplaneoffset(3);
        
%     fieldpattern = cell(3,1);
    Efieldpattern = cell(3,1);
    fieldpatternmatrixsize = cell(3,1);
    fieldpatternmatrixsize{1} = [fpms(1) fpms(2) 1];
    fieldpatternmatrixsize{2} = [fpms(1) 1 fpms(3)];
    fieldpatternmatrixsize{3} = [1 fpms(2) fpms(3)];
    Efieldpattern{1} = zeros([fieldpatternmatrixsize{1} 3 size(whichvoxels,1)]); % for each voxel there is matrix with field values
    Efieldpattern{2} = zeros([fieldpatternmatrixsize{2} 3 size(whichvoxels,1)]); % for each voxel there is matrix with field values
    Efieldpattern{3} = zeros([fieldpatternmatrixsize{3} 3 size(whichvoxels,1)]); % for each voxel there is matrix with field values
    Bfieldpattern = Efieldpattern;
    if s_opts.compute_coil_field_patterns
        Efieldpattern_coil = Efieldpattern;
        Bfieldpattern_coil = Efieldpattern;
    else
        Efieldpattern_coil = Efieldpattern;
        Bfieldpattern_coil = Efieldpattern;
    end
    mask_r1 = cell(3,1);
    mask_r2 = cell(3,1);
    mask_r3 = cell(3,1);
    mask_r4 = cell(3,1);

    % Field patterns are always calculated for all layers
    fieldpatternfov(1) = 2*g_opts.radius1; 
    fieldpatternfov(2) = 2*g_opts.radius1;
    fieldpatternfov(3) = 2*g_opts.radius1;
    
    fieldpatternfovstart_x = -fieldpatternfov(1)/2;
    fieldpatternfovend_x = fieldpatternfov(1)/2;
    fieldpatternfovstart_y = -fieldpatternfov(2)/2;
    fieldpatternfovend_y = fieldpatternfov(2)/2;
    fieldpatternfovstart_z = -fieldpatternfov(3)/2;
    fieldpatternfovend_z = fieldpatternfov(3)/2;
    if mod(fpms(1),2), 
        delta_fieldx = fieldpatternfov(1)/(fpms(1)-1);
    else
        delta_fieldx = fieldpatternfov(1)/fpms(1);
    end
    if mod(fpms(2),2), 
        delta_fieldy = fieldpatternfov(2)/(fpms(2)-1);
    else
        delta_fieldy = fieldpatternfov(2)/fpms(2);
    end
    if mod(fpms(3),2), 
        delta_fieldz = fieldpatternfov(3)/(fpms(3)-1);
    else
        delta_fieldz = fieldpatternfov(3)/fpms(3);
    end
    x_axis = fieldpatternfovstart_x:delta_fieldx:(fieldpatternfovend_x-~mod(fpms(1),2)*delta_fieldx);
    y_axis = fieldpatternfovstart_y:delta_fieldy:(fieldpatternfovend_y-~mod(fpms(2),2)*delta_fieldy);
    z_axis = fieldpatternfovstart_z:delta_fieldz:(fieldpatternfovend_z-~mod(fpms(3),2)*delta_fieldz);
    
    % make voxel distribution symmetric with respect to the origin if matrix size is even
    if mod(fpms(1),2) == 0
        x_axis = x_axis + ((fieldpatternfov(1)/fpms(1))/2);
    end
    if mod(fpms(2),2) == 0
        y_axis = y_axis + ((fieldpatternfov(2)/fpms(2))/2);
    end
    if mod(fpms(3),2) == 0
        z_axis = z_axis + ((fieldpatternfov(3)/fpms(3))/2);
    end
            
    if fpms(1) == 1,
        x_axis = 0;
    end
    if fpms(2) == 1,
        y_axis = 0;
    end
    if fpms(3) == 1,
        z_axis = 0;
    end
        fieldx = cell(3,1);
        fieldy = fieldx;
        fieldz = fieldx;
        
        fieldrho = fieldx;
        cosfieldtheta = fieldx;
        fieldtheta = fieldx; 
        fieldphi = fieldx;
        sinfieldtheta = fieldx;
        sinfieldphi = fieldx;
        cosfieldphi = fieldx;
        
        [fieldx{1},fieldy{1},fieldz{1}] = ndgrid(x_axis+fieldxoffset,y_axis+fieldyoffset,fieldzoffset);
        fieldrho{1} = sqrt(fieldx{1}.^2 + fieldy{1}.^2 + fieldz{1}.^2);
        cosfieldtheta{1} = fieldz{1}./fieldrho{1};
        cosfieldtheta{1}(isnan(cosfieldtheta{1})) = 0;        %** NaN at r=0  
        fieldtheta{1} = acos(cosfieldtheta{1});               % theta
        fieldphi{1} = atan2(fieldy{1},fieldx{1});              % phi
        sinfieldtheta{1} = sin(fieldtheta{1});              % sin[theta]
        sinfieldphi{1} = sin(fieldphi{1});                  % sin[phi]
        cosfieldphi{1} = cos(fieldphi{1});                  % cos[phi]
        cotfieldtheta{1} = cot(fieldtheta{1});              % cot[theta] ** Inf at theta=0
        cscfieldtheta{1} = cot(fieldtheta{1});              % csc[theta] ** Inf at theta=0
        cos2fieldtheta{1} = cosfieldtheta{1}.*cosfieldtheta{1};   % csc[theta] ** Inf at theta=0
        
        mask_r4{1} = fieldrho{1} < radius3; % REGION 4
        mask_r3{1} = (fieldrho{1} < radius2) & (fieldrho{1} >= radius3); % REGION 3
        mask_r2{1} = (fieldrho{1} < radius1) & (fieldrho{1} >= radius2); % REGION 2
        mask_r1{1} = ~(mask_r2{1} + mask_r3{1} + mask_r4{1}); % REGION 1
        if 0
            figure;
            subplot(2,2,1); imshow(squeeze(mask_r4{1}),[]);
            subplot(2,2,2); imshow(squeeze(mask_r3{1}),[]);
            subplot(2,2,3); imshow(squeeze(mask_r2{1}),[]);
            subplot(2,2,4); imshow(squeeze(mask_r1{1}),[]);
        end
        [fieldx{2},fieldy{2},fieldz{2}] = ndgrid(x_axis+fieldxoffset,fieldyoffset,z_axis+fieldzoffset);
%         % TO AVOID THETA == 0
%         if (fieldxoffset==0) && (fieldyoffset==0)
%             [fieldx{2},fieldy{2},fieldz{2}] = ndgrid(x_axis+fieldxoffset+0.0015,fieldyoffset+0.0015,z_axis+fieldzoffset);
%         elseif (fieldxoffset==0) && ~(fieldyoffset==0)
%             [fieldx{2},fieldy{2},fieldz{2}] = ndgrid(x_axis+fieldxoffset+0.0015,fieldyoffset,z_axis+fieldzoffset);
%         elseif ~(fieldxoffset==0) && (fieldyoffset==0)
%             [fieldx{2},fieldy{2},fieldz{2}] = ndgrid(x_axis+fieldxoffset,fieldyoffset+0.0015,z_axis+fieldzoffset);
%         else
%             [fieldx{2},fieldy{2},fieldz{2}] = ndgrid(x_axis+fieldxoffset,fieldyoffset,z_axis+fieldzoffset);
%         end
        fieldrho{2} = sqrt(fieldx{2}.^2 + fieldy{2}.^2 + fieldz{2}.^2);
        cosfieldtheta{2} = fieldz{2}./fieldrho{2};
        cosfieldtheta{2}(isnan(cosfieldtheta{2})) = 0;        %** NaN at r=0  
        fieldtheta{2} = acos(cosfieldtheta{2});               % theta
        fieldphi{2} = atan2(fieldy{2},fieldx{2});              % phi
        sinfieldtheta{2} = sin(fieldtheta{2});              % sin[theta]
        sinfieldphi{2} = sin(fieldphi{2});                  % sin[phi]
        cosfieldphi{2} = cos(fieldphi{2});                  % cos[phi]
        cotfieldtheta{2} = cot(fieldtheta{2});              % cot[theta] ** Inf at theta=0
        cscfieldtheta{2} = cot(fieldtheta{2});              % csc[theta] ** Inf at theta=0
        cos2fieldtheta{2} = cosfieldtheta{2}.*cosfieldtheta{2};   % csc[theta] ** Inf at theta=0
        
        mask_r4{2} = fieldrho{2} < radius3; % REGION 4
        mask_r3{2} = (fieldrho{2} < radius2) & (fieldrho{2} >= radius3); % REGION 3
        mask_r2{2} = (fieldrho{2} < radius1) & (fieldrho{2} >= radius2); % REGION 2
        mask_r1{2} = ~(mask_r2{2} + mask_r3{2} + mask_r4{2}); % REGION 1
        if 0
            figure;
            subplot(2,2,1); imshow(squeeze(mask_r4{2}),[]);
            subplot(2,2,2); imshow(squeeze(mask_r3{2}),[]);
            subplot(2,2,3); imshow(squeeze(mask_r2{2}),[]);
            subplot(2,2,4); imshow(squeeze(mask_r1{2}),[]);
        end
        [fieldx{3},fieldy{3},fieldz{3}] = ndgrid(fieldxoffset,y_axis+fieldyoffset,z_axis+fieldzoffset);
%         % TO AVOID THETA == 0
%         if (fieldxoffset==0) && (fieldyoffset==0)
%             [fieldx{3},fieldy{3},fieldz{3}] = ndgrid(fieldxoffset+0.0015,y_axis+fieldyoffset+0.0015,z_axis+fieldzoffset);
%         elseif (fieldxoffset==0) && ~(fieldyoffset==0)
%             [fieldx{3},fieldy{3},fieldz{3}] = ndgrid(fieldxoffset+0.0015,y_axis+fieldyoffset,z_axis+fieldzoffset);
%         elseif ~(fieldxoffset==0) && (fieldyoffset==0)
%             [fieldx{3},fieldy{3},fieldz{3}] = ndgrid(fieldxoffset,y_axis+fieldyoffset+0.0015,z_axis+fieldzoffset);
%         else
%             [fieldx{3},fieldy{3},fieldz{3}] = ndgrid(fieldxoffset,y_axis+fieldyoffset,z_axis+fieldzoffset);
%         end        
        
        fieldrho{3} = sqrt(fieldx{3}.^2 + fieldy{3}.^2 + fieldz{3}.^2);
        cosfieldtheta{3} = fieldz{3}./fieldrho{3};
        cosfieldtheta{3}(isnan(cosfieldtheta{3})) = 0;        %** NaN at r=0  
        fieldtheta{3} = acos(cosfieldtheta{3});               % theta
        fieldphi{3} = atan2(fieldy{3},fieldx{3});              % phi
        sinfieldtheta{3} = sin(fieldtheta{3});              % sin[theta]
        sinfieldphi{3} = sin(fieldphi{3});                  % sin[phi]
        cosfieldphi{3} = cos(fieldphi{3});                  % cos[phi]
        cotfieldtheta{3} = cot(fieldtheta{3});              % cot[theta] ** Inf at theta=0
        cscfieldtheta{3} = cot(fieldtheta{3});              % csc[theta] ** Inf at theta=0
        cos2fieldtheta{3} = cosfieldtheta{3}.*cosfieldtheta{3};   % csc[theta] ** Inf at theta=0
        
%         [smtheta_1, smtheta_2, smtheta_3] = find(abs(fieldtheta{3}) < 0.01);
        
        mask_r4{3} = fieldrho{3} < radius3; % REGION 4
        mask_r3{3} = (fieldrho{3} < radius2) & (fieldrho{3} >= radius3); % REGION 3
        mask_r2{3} = (fieldrho{3} < radius1) & (fieldrho{3} >= radius2); % REGION 2
        mask_r1{3} = ~(mask_r2{3} + mask_r3{3} + mask_r4{3}); % REGION 1
        if 0
            figure;
            subplot(2,2,1); imshow(squeeze(mask_r4{3}),[]);
            subplot(2,2,2); imshow(squeeze(mask_r3{3}),[]);
            subplot(2,2,3); imshow(squeeze(mask_r2{3}),[]);
            subplot(2,2,4); imshow(squeeze(mask_r1{3}),[]);
        end
        
else
    Efieldpattern = [];
    Bfieldpattern = [];
    Efieldpattern_coil = [];
    Bfieldpattern_coil = [];
end

% -- initialize SNR and g-factor -- 

snr_indx = find(snr_fov == 1);
[snr_row,snr_col] = find(snr_fov == 1);
mask_snr = mask_snr(min(snr_row):max(snr_row),min(snr_col):max(snr_col));
% NOTE: mask_snr and snr_fov have the same number of nonzero elements, but the snr_fov has the dimensions of the full matrix_size (i.e., including all regions)
tempindex_r = max(snr_row)-(nf/2);
tempindex_c = max(snr_col)-(np/2);
snr_ult = zeros(2*tempindex_r,2*tempindex_c);
g_ult = snr_ult;
if s_opts.compute_coilSNR_flag
    snr_coil = snr_ult;
    g_coil = g_ult;
else
    snr_coil = [];
    g_coil = [];
    psi_coil = [];
end

% -- initialize EM fields -- 
if s_opts.save_emfields_ult_r4_flag || s_opts.save_emfields_coil_r4_flag || s_opts.save_emfields_ult_r3_flag || s_opts.save_emfields_coil_r3_flag || s_opts.save_emfields_ult_r2_flag || s_opts.save_emfields_coil_r2_flag 
    save_emfields_flag = 1;
    if s_opts.save_emfields_ult_r4_flag || s_opts.save_emfields_ult_r3_flag || s_opts.save_emfields_ult_r2_flag
        efield_ult_set = zeros(numbasis,2*tempindex_r,2*tempindex_c,3);
        bfield_ult_set = zeros(numbasis,2*tempindex_r,2*tempindex_c,3);
    else
        efield_ult_set_fullfov = [];
        bfield_ult_set_fullfov = [];
    end
    if s_opts.compute_coilSNR_flag && (s_opts.save_emfields_coil_r4_flag || s_opts.save_emfields_coil_r3_flag || s_opts.save_emfields_coil_r2_flag)
        efield_coil_set = zeros(ncoils,2*tempindex_r,2*tempindex_c,3);
        bfield_coil_set = zeros(ncoils,2*tempindex_r,2*tempindex_c,3);
    else
        efield_coil_set_fullfov = [];
        bfield_coil_set_fullfov = [];
    end
else
    save_emfields_flag = 0;
    efield_ult_set_fullfov = [];
    efield_coil_set_fullfov = [];
    bfield_ult_set_fullfov = [];
    bfield_coil_set_fullfov = [];
end

%--------------------------------------------
%        Ultimate intrinsic SNR scaling
%--------------------------------------------

Nproton = 6.691e28;     % number of protons per unit volume in water
gyromag = 2.68e8;       % gyromagnetic ratio for protons [rad/T/sec]
h = 6.626e-34;          % Planck's constant [J sec]
hbar = h/2/pi;
Ispin = 1/2;
k_B = 1.3806503e-23;    % Boltzmann's constant
T = 310;
M_0 = Nproton*((gyromag*hbar)^2)*Ispin*(Ispin+1)*fieldstrength/(3*k_B*T);

usnr_num = omega*M_0/sqrt(8*k_B*T);

%** experimental SNR scaling is defined in the batch file 
snr_num = usnr_num*expsnr_num;

% ---------------------------------- %
%         Simulation begins
% ---------------------------------- %
%
% ORGANIZATION:
%        - calculate SNR, EM field at voxels for which the SNR is calculated
%        - calculate current patterns at specified voxel/s (among those for which SNR is calculated)
%        - calculate field patterns using weights for optimal SNR at specific voxel locations
%        - calculate EM field for regions 2 and 3 if not included in the SNR calculation
%        - EM field and radiation losses for region 1
%-----------------------------------

% -- constant factors depending on radius -- %
% k_0_rad_a = k_0*radius3;
k_0_rad_b = k_0*currentradius;
% besselnorm_rad_a = sqrt(pi/(2*k_0_rad_a));
besselnorm_rad_b = sqrt(pi/(2*k_0_rad_b));

ele_scaling = omega*mu*k_0*currentradius*currentradius; % scaling factor for the Electric field (RL fixed on 6/11/2012)

mag_scaling_r4 = -1i*mu*k_4*k_0*currentradius*currentradius; % scaling factor for the Magnetic field (RL fixed on 6/11/2012)
mag_scaling_r3 = -1i*mu*k_3*k_0*currentradius*currentradius; 
mag_scaling_r2 = -1i*mu*k_2*k_0*currentradius*currentradius; 
mag_scaling_r1 = -1i*mu*k_0*k_0*currentradius*currentradius;

% -- loop through voxel positions -- 
matsize_snr = size(snr_ult); % Need to use a matrix fitted around SNR fov for the accelerated case
ind_skip = floor(matsize_snr./acceleration);

disp('  looping through voxel positions...')
voxel_count = 0;
for ind_1 = 1:ind_skip(1), 
    ind_1_set = ind_1:ind_skip(1):matsize_snr(1);
    for ind_2 = 1:ind_skip(2),
        disp(['    voxel (' num2str(ind_1) ',' num2str(ind_2) ')'])
        ind_2_set = ind_2:ind_skip(2):matsize_snr(2);
             
        voxel_count = voxel_count + 1; % used to avoid calculating noise power more than once
                  
        xset = x_fov(ind_1_set+(min(snr_row)-1),ind_2_set+(min(snr_col)-1));
        xset = xset(:).';
        yset = y_fov(ind_1_set+(min(snr_row)-1),ind_2_set+(min(snr_col)-1));
        yset = yset(:).';
        zset = z_fov(ind_1_set+(min(snr_row)-1),ind_2_set+(min(snr_col)-1));
        zset = zset(:).';

        disp(['    *** x,y,z [cm] = (' num2str(100*xset) ','  num2str(100*yset) ',' num2str(100*zset) ')'])

        % Convert in spherical coordinates
        rset = sqrt(xset.^2 + yset.^2 + zset.^2);   % rho
        costhetaset = zset./rset;                   % cos[theta]
        costhetaset(isnan(costhetaset)) = 0;        %** NaN at r=0  
        thetaset = acos(costhetaset);               % theta
        phiset = atan2(yset,xset);                  % phi

        sinthetaset = sin(thetaset);                % sin[theta]
        sinphiset = sin(phiset);                    % sin[phi]
        cosphiset = cos(phiset);                    % cos[phi]

        cotthetaset = cot(thetaset);                % cot[theta] ** Inf at theta=0
        cscthetaset = csc(thetaset);                % csc[theta] ** Inf at theta=0
        cos2thetaset = costhetaset.*costhetaset;    % cos[theta]^2
        
        % unit vectors conversion
        rhat_x = sinthetaset.*cosphiset;
        rhat_y = sinthetaset.*sinphiset;
        rhat_z = costhetaset;
        
        if snr_radius == 3, % ONLY INNERMOST LAYER
            krset = k_4*rset;
        else
            krset = zeros(size(rset));
            for irho = 1:length(rset)
                if rset(irho) < radius3
                    krset(irho) = k_4*rset(irho);
                    disp('radius < R3')
                elseif (rset(irho) < radius2) && (rset(irho) >= radius3)
                    krset(irho) = k_3*rset(irho);
                    disp('R3 < radius < R2')
                elseif (rset(irho) < radius1) && (rset(irho) >= radius2)
                    krset(irho) = k_2*rset(irho);
                    disp('R2 < radius < R1')
                else
                    krset(irho) = k_0*rset(irho);
                    disp('radius > R1')
                end
            end
        end
        
        besselnorm = sqrt(pi./(2.*krset));
        % ----------------------------------------- %
        
        % generate sensitivity matrix for each mode
        counter = 1;
        counter_modes = 1;
        X_full = zeros(numbasis,length(ind_1_set)*length(ind_2_set));
        X_P_full = X_full.';
        
        if s_opts.compute_coilSNR_flag,
            % Initialize quantities needed for circular coil computations
            X_modes = zeros(((lmax + 1)^2 - 1),length(ind_1_set)*length(ind_2_set));
            P_modes = zeros(((lmax + 1)^2 - 1),1);
%             P_L_modes_r4 = zeros(((lmax + 1)^2 - 1),1);
%             P_L_modes_r3 = zeros(((lmax + 1)^2 - 1),1);
%             P_L_modes_r2 = zeros(((lmax + 1)^2 - 1),1);
            P_L_modes = zeros(((lmax + 1)^2 - 1),1);
            P_R_modes = zeros(((lmax + 1)^2 - 1),1);
        end
        
        if save_emfields_flag,
            all_E_phase = zeros(numbasis,length(ind_1_set)*length(ind_2_set),3);
            all_B_phase = zeros(numbasis,length(ind_1_set)*length(ind_2_set),3);
        end
        if s_opts.compute_coilSNR_flag && save_emfields_flag,
            all_E_phase_modes = zeros(((lmax + 1)^2 - 1),length(ind_1_set)*length(ind_2_set),3);
            all_E_phase_coil = zeros(ncoils,length(ind_1_set)*length(ind_2_set),3);
            all_B_phase_modes = zeros(((lmax + 1)^2 - 1),length(ind_1_set)*length(ind_2_set),3);
            all_B_phase_coil = zeros(ncoils,length(ind_1_set)*length(ind_2_set),3);
        end
        
        T_matrix_set = cell(numbasis);
        
        for l = 1:lmax, % for l=0 the T matrix is empty and there would be no contributions 
            %             disp(['      l = ' num2str(l)])
            lnorm = sqrt(l*(l+1));
            legendrenorm = sqrt((2*l + 1)/(4*pi));
            legendrenorm_minus1 = sqrt((2*l - 1)/(4*pi));
            legendrefunctions = legendre(l,costhetaset,'sch'); % # row is m = 0,...l ; # col is R
            if l>0,
                legendrefunctions_lminus1 = legendre(l-1,costhetaset,'sch');
                legendrefunctions_lminus1 = [legendrefunctions_lminus1; zeros(size(costhetaset))];
            end

            % -- these two are for rho = sphere radius --
%             H_l_1_k_0_rad_a = besselnorm_rad_a*besselh(l+0.5,1,k_0_rad_a);
%             H_l_1_k_0_rad_a_prime = -k_0_rad_a*besselnorm_rad_a*besselh(l+1+0.5,1,k_0_rad_a) + (l + 1)*H_l_1_k_0_rad_a;
            % -----------------------------
          
            H_l_1_k_0_rad_b = besselnorm_rad_b*besselh(l+0.5,1,k_0_rad_b);
            H_l_1_k_0_rad_b_prime = -k_0_rad_b*besselnorm_rad_b*besselh(l+1+0.5,1,k_0_rad_b) + (l + 1)*H_l_1_k_0_rad_b;
            
            J_l_kr = besselnorm.*besselj(l+0.5,krset);
            J_l_kr(rset<eps) = (l==0);                            %** only j_0 survives at r=0
            J_lplus1_kr = besselnorm.*besselj(l+1+0.5,krset);
            J_lplus1_kr(rset<eps) = 0;                           %** j_l+1 vanishes at r=0 for l>=0
            J_l_kr_prime = -krset.*J_lplus1_kr + (l + 1)*J_l_kr;
            
            H_l_1_kr = besselnorm.*besselh(l+0.5,1,krset);
            H_l_1_kr_prime = -krset.*besselnorm.*besselh(l+1+0.5,1,krset) + (l + 1)*H_l_1_kr;
            
            [R_H_P_1, R_H_F_1, R_V_P_1, R_V_F_1, T_H_P_1, T_H_F_1, T_V_P_1, T_V_F_1] = ...
                compute_reflection_and_transmission_coef(l,k_0,k_2,g_opts.radius1);
            
            [R_H_P_2, R_H_F_2, R_V_P_2, R_V_F_2, T_H_P_2, T_H_F_2, T_V_P_2, T_V_F_2] = ...
                compute_reflection_and_transmission_coef(l,k_2,k_3,g_opts.radius2);
            
            [R_H_P_3, R_H_F_3, R_V_P_3, R_V_F_3, T_H_P_3, T_H_F_3, T_V_P_3, T_V_F_3] = ...
                compute_reflection_and_transmission_coef(l,k_3,k_4,radius3);
            
            B_m_11 = -( T_H_P_2*(T_H_P_1*R_H_F_1 + T_H_F_1*R_H_F_2) + R_H_F_3*T_H_F_2*(T_H_F_1 + T_H_P_1*R_H_P_2*R_H_F_1) )/...
                ( T_H_P_2*(T_H_P_1 + T_H_F_1*R_H_F_2*R_H_P_1) + R_H_F_3*T_H_F_2*(T_H_P_1*R_H_P_2 + T_H_F_1*R_H_P_1) );
            
            
            B_n_11 = -( T_V_P_2*(T_V_P_1*R_V_F_1 + T_V_F_1*R_V_F_2) + R_V_F_3*T_V_F_2*(T_V_F_1 + T_V_P_1*R_V_P_2*R_V_F_1) )/...
                ( T_V_P_2*(T_V_P_1 + T_V_F_1*R_V_F_2*R_V_P_1) + R_V_F_3*T_V_F_2*(T_V_P_1*R_V_P_2 + T_V_F_1*R_V_P_1) );
            
            B_m_21 = (1/T_H_F_1)*(B_m_11 + R_H_F_1);
            
            B_n_21 = (1/T_V_F_1)*(B_n_11 + R_V_F_1);
            
            D_m_21 = (1/T_H_P_1)*(1+ R_H_P_1*B_m_11);
            
            D_n_21 = (1/T_V_P_1)*(1+ R_V_P_1*B_n_11);
            
            B_m_31 = (1/T_H_F_2)*(B_m_21 + R_H_F_2*D_m_21);
            
            B_n_31 = (1/T_V_F_2)*(B_n_21 + R_V_F_2*D_n_21);
            
            D_m_31 = (1/T_H_P_2)*(R_H_P_2*B_m_21 + D_m_21);
            
            D_n_31 = (1/T_V_P_2)*(R_V_P_2*B_n_21 + D_n_21);
            
            D_m_41 = (1/T_H_P_3)*(R_H_P_3*B_m_31 + D_m_31);
            
            D_n_41 = (1/T_V_P_3)*(R_V_P_3*B_n_31 + D_n_31);
            
            
            % T matrix does not depend on the expansion parameter m
            Ttemp = [H_l_1_k_0_rad_b                     0
                            0            (1/k_0_rad_b)*H_l_1_k_0_rad_b_prime].';
            
            for m = -l:l,
                %                 disp(['        m = ' num2str(m)])
%                 lmul = sqrt( (2*l + 1)*(l^2 - m^2)/(4*pi) );

                T = ((-1)^(1-m))*Ttemp;
                
%             disp(['l = ' num2str(l)])
%             disp(['    m = ' num2str(m)])
%             disp(['        factor = ' num2str((-1)^(1-m))])
%             disp(['T11 = ' num2str(T(1,1)) ', T22 = ' num2str(T(2,2))])
                

                lmul = sqrt((2*l + 1)*(l^2 - m^2)/(2*l - 1));
                
                if m > 0
                    Y_l_m = ((-1)^m)*legendrenorm*(1/sqrt(2))*legendrefunctions((m+1),:).*exp(1i*m*phiset);
                elseif m == 0
                    Y_l_m = legendrenorm*legendrefunctions(1,:);
                else % using Y_(l,-m) = (-1)^m*conj[Y_(l,m)]
                    Y_l_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm*(1/sqrt(2))*legendrefunctions((abs(m)+1),:).*exp(1i*abs(m)*phiset));
                end
                S_E = zeros(size(Y_l_m));
                S_M = S_E;
                if m > 0
                    Y_lminus1_m = ((-1)^m)*legendrenorm_minus1*(1/sqrt(2))*legendrefunctions_lminus1(m+1,:).*exp(1i*m*phiset);
                elseif m == 0
                    Y_lminus1_m = legendrenorm_minus1*legendrefunctions_lminus1(1,:);
                else
                    Y_lminus1_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm_minus1*(1/sqrt(2))*legendrefunctions_lminus1((abs(m)+1),:).*exp(1i*abs(m)*phiset));
                end
                   
                X_x = (1/lnorm).*((-m*cosphiset+1i*l*sinphiset).*cotthetaset.*Y_l_m ...
                    - 1i*lmul*cscthetaset.*sinphiset.*Y_lminus1_m);
                X_y = (1/lnorm).*((-m*sinphiset-1i*l*cosphiset).*cotthetaset.*Y_l_m ...
                    + 1i*lmul*cscthetaset.*cosphiset.*Y_lminus1_m);
                X_z = (1/lnorm).*m.*Y_l_m;

                r_cross_X_x = (1/lnorm).*((m*sinphiset+1i*l*cos2thetaset.*cosphiset).*cscthetaset.*Y_l_m ...
                    - 1i*lmul*cotthetaset.*cosphiset.*Y_lminus1_m);
                r_cross_X_y = (1/lnorm).*((-m*cosphiset+1i*l*cos2thetaset.*sinphiset).*cscthetaset.*Y_l_m ...
                    - 1i*lmul*cotthetaset.*sinphiset.*Y_lminus1_m);
                r_cross_X_z = (-1i/lnorm).*(l*costhetaset.*Y_l_m ...
                    - lmul.*Y_lminus1_m);
                
                if rset < eps
                    M_x = 0;
                    M_y = 0;
                    M_z = 0;
                    if l == 1,
                        switch m
                            case -1
                                N_x=1i*sqrt(1/(6*pi))*(1/sqrt(2));
                                N_y=1i*sqrt(1/(6*pi))*(-1i/sqrt(2));
                                N_z=0;
                            case 0
                                N_x=0;
                                N_y=0;
                                N_z=1i*sqrt(1/(6*pi));
                            case 1
                                N_x=1i*sqrt(1/(6*pi))*(-1/sqrt(2));
                                N_y=1i*sqrt(1/(6*pi))*(-1i/sqrt(2));
                                N_z=0;
                        end
                    else
                        N_x=0;
                        N_y=0;
                        N_z=0;
                    end
                else
                    M_x = J_l_kr.*X_x;
                    M_y = J_l_kr.*X_y;
                    M_z = J_l_kr.*X_z;
                    
                    N_x = (1./krset).*J_l_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_x;
                    N_y = (1./krset).*J_l_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_y;
                    N_z = (1./krset).*J_l_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_z;
                end
                
                M3_x = H_l_1_kr.*X_x;
                M3_y = H_l_1_kr.*X_y;
                M3_z = H_l_1_kr.*X_z;
                
                N3_x = (1./krset).*H_l_1_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_x;
                N3_y = (1./krset).*H_l_1_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_y;
                N3_z = (1./krset).*H_l_1_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_z;

                % create sensitivity matrices using B- receive fields
                
                S_M = zeros(1,length(rset));
                S_E = zeros(1,length(rset));
                for irho = 1:length(rset)
                    if rset(irho) < radius3
                        S_M(irho) = mag_scaling_r4*D_m_41*(N_x(irho) - 1i*N_y(irho));
                        S_E(irho) = mag_scaling_r4*D_n_41*(M_x(irho) - 1i*M_y(irho));
                    elseif (rset(irho) < radius2) && (rset(irho) >= radius3)
                        S_M(irho) = mag_scaling_r3*(B_m_31*(N3_x(irho) - 1i*N3_y(irho)) + D_m_31*(N_x(irho) - 1i*N_y(irho)));
                        S_E(irho) = mag_scaling_r3*(B_n_31*(M3_x(irho) - 1i*M3_y(irho)) + D_n_31*(M_x(irho) - 1i*M_y(irho)));
                    elseif (rset(irho) < radius1) && (rset(irho) >= radius2)
                        S_M(irho) = mag_scaling_r2*(B_m_21*(N3_x(irho) - 1i*N3_y(irho)) + D_m_21*(N_x(irho) - 1i*N_y(irho)));
                        S_E(irho) = mag_scaling_r2*(B_n_21*(M3_x(irho) - 1i*M3_y(irho)) + D_n_21*(M_x(irho) - 1i*M_y(irho)));
                    else
                        S_M(irho) = 0;
                        S_E(irho) = 0;
                    end
                end
                
                S = [S_M ; S_E];

                T_hat = T;
                if whichcurrents == 1,
                    T_hat = [1 0]*T;
                end
                if whichcurrents == 2,
                    T_hat = [0 1]*T;
                end
                
                T_matrix_set{counter_modes} = T_hat;
                
                % ------------ %
                % form X = T*S
                X = T_hat*S;

                % -- SAVE E AND B FIELDS BASED ON REGION -- %
                if save_emfields_flag
                    if snr_radius == 3, % then all aliasing voxels belong to the same region
                        tempmat = [D_m_41; D_n_41];
                        tempmat = repmat(tempmat,1,length(M_x));
                        % Create EM field components
                        E_x = ele_scaling*(T_hat*(tempmat.*[M_x; N_x]));
                        E_y = ele_scaling*(T_hat*(tempmat.*[M_y; N_y]));
                        E_z = ele_scaling*(T_hat*(tempmat.*[M_z; N_z]));
                        
                        B_x = mag_scaling_r4*(T_hat*(tempmat.*[N_x; M_x]));
                        B_y = mag_scaling_r4*(T_hat*(tempmat.*[N_y; M_y]));
                        B_z = mag_scaling_r4*(T_hat*(tempmat.*[N_z; M_z]));
                        clear tempmat
                    else 
%                         E_x = zeros(2,length(M_x));
%                         E_y = zeros(2,length(M_x));
%                         E_z = zeros(2,length(M_x));
%                         B_x = zeros(2,length(M_x));
%                         B_y = zeros(2,length(M_x));
%                         B_z = zeros(2,length(M_x));
                        E_x = zeros((counterincrement+1),length(M_x)); % RL 30 Jan 2015
                        E_y = zeros((counterincrement+1),length(M_x));
                        E_z = zeros((counterincrement+1),length(M_x));
                        B_x = zeros((counterincrement+1),length(M_x));
                        B_y = zeros((counterincrement+1),length(M_x));
                        B_z = zeros((counterincrement+1),length(M_x));
                        for irho = 1:length(rset) % if doing parallel imaging, the aliasing voxels may belong to different regions
                            if rset(irho) < radius3 % REGION 4
                                E_x(:,irho) = ele_scaling*(T_hat*([D_m_41; D_n_41].*[M_x(irho); N_x(irho)]));
                                E_y(:,irho) = ele_scaling*(T_hat*([D_m_41; D_n_41].*[M_y(irho); N_y(irho)]));
                                E_z(:,irho) = ele_scaling*(T_hat*([D_m_41; D_n_41].*[M_z(irho); N_z(irho)]));
                                B_x(:,irho) = mag_scaling_r4*(T_hat*([D_m_41; D_n_41].*[N_x(irho); M_x(irho)]));
                                B_y(:,irho) = mag_scaling_r4*(T_hat*([D_m_41; D_n_41].*[N_y(irho); M_y(irho)]));
                                B_z(:,irho) = mag_scaling_r4*(T_hat*([D_m_41; D_n_41].*[N_z(irho); M_z(irho)]));
                            elseif (rset(irho) < radius2) && (rset(irho) >= radius3) % REGION 3
                                E_x(:,irho) = ele_scaling*(T_hat*[ (B_m_31*M3_x(irho) + D_m_31*M_x(irho)); (B_n_31*N3_x(irho) + D_n_31*N_x(irho)) ]);
                                E_y(:,irho) = ele_scaling*(T_hat*[ (B_m_31*M3_y(irho) + D_m_31*M_y(irho)); (B_n_31*N3_y(irho) + D_n_31*N_y(irho)) ]);
                                E_z(:,irho) = ele_scaling*(T_hat*[ (B_m_31*M3_z(irho) + D_m_31*M_z(irho)); (B_n_31*N3_z(irho) + D_n_31*N_z(irho)) ]);
                                B_x(:,irho) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_x(irho) + D_m_31*N_x(irho)); (B_n_31*M3_x(irho) + D_n_31*M_x(irho)) ]);
                                B_y(:,irho) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_y(irho) + D_m_31*N_y(irho)); (B_n_31*M3_y(irho) + D_n_31*M_y(irho)) ]);
                                B_z(:,irho) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_z(irho) + D_m_31*N_z(irho)); (B_n_31*M3_z(irho) + D_n_31*M_z(irho)) ]);
                            elseif (rset(irho) < radius1) && (rset(irho) >= radius2) % REGION 2
                                E_x(:,irho) = ele_scaling*(T_hat*[ (B_m_21*M3_x(irho) + D_m_21*M_x(irho)); (B_n_21*N3_x(irho) + D_n_21*N_x(irho)) ]);
                                E_y(:,irho) = ele_scaling*(T_hat*[ (B_m_21*M3_y(irho) + D_m_21*M_y(irho)); (B_n_21*N3_y(irho) + D_n_21*N_y(irho)) ]);
                                E_z(:,irho) = ele_scaling*(T_hat*[ (B_m_21*M3_z(irho) + D_m_21*M_z(irho)); (B_n_21*N3_z(irho) + D_n_21*N_z(irho)) ]);
                                B_x(:,irho) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_x(irho) + D_m_21*N_x(irho)); (B_n_21*M3_x(irho) + D_n_21*M_x(irho)) ]);
                                B_y(:,irho) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_y(irho) + D_m_21*N_y(irho)); (B_n_21*M3_y(irho) + D_n_21*M_y(irho)) ]);
                                B_z(:,irho) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_z(irho) + D_m_21*N_z(irho)); (B_n_21*M3_z(irho) + D_n_21*M_z(irho)) ]);
                            else % REGION 1 (outside the FOV, calculated separately from SNR) 
%                                 E_x(:,irho) = [0; 0];
%                                 E_y(:,irho) = [0; 0];
%                                 E_z(:,irho) = [0; 0];
%                                 B_x(:,irho) = [0; 0];
%                                 B_y(:,irho) = [0; 0];
%                                 B_z(:,irho) = [0; 0];
                                E_x(:,irho) = zeros((counterincrement+1),1); % RL 30 Jan 2015
                                E_y(:,irho) = zeros((counterincrement+1),1);
                                E_z(:,irho) = zeros((counterincrement+1),1);
                                B_x(:,irho) = zeros((counterincrement+1),1);
                                B_y(:,irho) = zeros((counterincrement+1),1);
                                B_z(:,irho) = zeros((counterincrement+1),1);
                            end
                        end
                    end
                    all_E_phase(counter:counter+counterincrement,:,:) = reshape([E_x E_y E_z],[(counterincrement+1) size(xset,2) 3]);
                    all_B_phase(counter:counter+counterincrement,:,:) = reshape([B_x B_y B_z],[(counterincrement+1) size(xset,2) 3]);
                end
                
                if voxel_count == 1 % Calculate noise covariance matrix only once
                    % ------------------------------------------------------------- %
                    % form noise correlation submatrix P (2 x 2 array for each l,m)
                    
                    % -- body losses component P_L
                    %                 P_L = (sigma/2)*(abs(omega*mu/k_4)^2)*...
                    %                     [ computePsiM(l, k_4, radius3)              0
                    %                     0               computePsiE(l, k_4, radius3)];
                    P_L_scaling = abs(omega*mu*k_0*currentradius*currentradius)^2;
                    if sigma_r4 == 0
                        P_L_r4 = [0 0; 0 0];
                    else
                        P_L_r4 = (sigma_r4/2)*P_L_scaling*...
                            [ (abs(D_m_41)^2)*computePsiM_allregions(l, k_4, 0, radius3, 1)              0
                            0                   (abs(D_n_41)^2)*computePsiE_allregions(l, k_4, 0, radius3, 1)];
                    end
                    if sigma_r3 == 0
                        P_L_r3 = [0 0; 0 0];
                    else
                        term1 = (abs(B_m_31)^2)*computePsiM_allregions(l, k_3, radius3 ,radius2, 2) + (abs(D_m_31)^2)*computePsiM_allregions(l, k_3, radius3 ,radius2, 1) + ...
                            (B_m_31*conj(D_m_31))*computePsiM_allregions(l, k_3, radius3 ,radius2, 4) + (D_m_31*conj(B_m_31))*computePsiM_allregions(l, k_3, radius3 ,radius2, 4);
                        
                        term2 = (abs(B_n_31)^2)*computePsiE_allregions(l, k_3, radius3 ,radius2, 2) + (abs(D_n_31)^2)*computePsiE_allregions(l, k_3, radius3 ,radius2, 1) + ...
                            (B_n_31*conj(D_n_31))*computePsiE_allregions(l, k_3, radius3 ,radius2, 4) + (D_n_31*conj(B_n_31))*computePsiE_allregions(l, k_3, radius3 ,radius2, 4);
                        
                        P_L_r3 = (sigma_r3/2)*P_L_scaling*...
                            [ term1             0
                            0             term2];
                        clear term1 term2
                    end
                    if sigma_r2 == 0
                        P_L_r2 = [0 0; 0 0];
                    else
                        term1 = (abs(B_m_21)^2)*computePsiM_allregions(l, k_2, radius2 ,radius1, 2) + (abs(D_m_21)^2)*computePsiM_allregions(l, k_2, radius2 ,radius1, 1) + ...
                            (B_m_21*conj(D_m_21))*computePsiM_allregions(l, k_2, radius2 ,radius1, 4) + (D_m_21*conj(B_m_21))*computePsiM_allregions(l, k_2, radius2 ,radius1, 4);
                        
                        term2 = (abs(B_n_21)^2)*computePsiE_allregions(l, k_2, radius2 ,radius1, 2) + (abs(D_n_21)^2)*computePsiE_allregions(l, k_2, radius2 ,radius1, 1) + ...
                            (B_n_21*conj(D_n_21))*computePsiE_allregions(l, k_2, radius2 ,radius1, 4) + (D_n_21*conj(B_n_21))*computePsiE_allregions(l, k_2, radius2 ,radius1, 4);
                        
                        P_L_r2 = (sigma_r2/2)*P_L_scaling*...
                            [ term1             0
                            0             term2];
                        clear term1 term2
                    end
                    P_L = P_L_r4 + P_L_r3 + P_L_r2;
                    
                    P = T_hat*P_L*T_hat';
                    
                    if s_opts.save_noise_contributions
                        P_L_set_r4(counter:counter+counterincrement) = diag(T_hat*P_L_r4*T_hat');
                        P_L_set_r3(counter:counter+counterincrement) = diag(T_hat*P_L_r3*T_hat');
                        P_L_set_r2(counter:counter+counterincrement) = diag(T_hat*P_L_r2*T_hat');
                    end
                    
                    % -- conductor losses component P_A
                    if s_opts.include_ult_conductor_losses || s_opts.include_coil_conductor_losses
                        P_A = (currentradius^2)/(2*sigma_coil*d_coil)*eye(2);
                        P_A_hat = P_A;
                        if whichcurrents == 1,
                            P_A_hat = P_A(1,1);
                        end
                        if whichcurrents == 2,
                            P_A_hat = P_A(2,2);
                        end
                        if s_opts.include_ult_conductor_losses
                            P = P + P_A_hat;
                            if s_opts.save_noise_contributions
                                P_A_set(counter:counter+counterincrement) = diag(P_A_hat);
                            end
                        end
                    end
                    
                    % -- @@@@@@@@@@@@@@@@@@ TO BE DONE!!! @@@@@@@@@@@@@@@@@@@@@@
                    % -- radiation losses component P_R % 
                    if s_opts.include_ult_radiation_losses || s_opts.include_coil_radiation_losses
                        %                     if m == 0
                        %                         P_R = sphere_calc_poynting(l,m,s_opts,g_opts,path_opts,plot_poynting_sphere_flag,whichcurrents,k_0,k_4,ele_scaling_poynting,mag_scaling_poynting);
                        %
                        %                     else
                        %                         P_R = 0;
                        %                     end
                        %                     H_l_1_k_0_rad_outer = besselnorm_rad_outer*besselh(l+0.5,1,k_0_rad_outer);
                        %
                        %                     n_l = (mu*2*pi)*sqrt((2*l+1)/(4*pi*l(l+1))*k_0^2*radius3^2*
                        %                     P_R_1=sqrt(1/mu*epsilon_0)/(2*mu*k_0^2)*n_l; % See Keltner 1991
                        
                        P_R = sphere_calc_poynting(l,m,s_opts,g_opts,path_opts,plot_poynting_sphere_flag,whichcurrents,k_0,k_4,ele_scaling_poynting,mag_scaling_poynting);
                        
                        P_R_hat = P_R;
                        if whichcurrents == 1,
                            P_R_hat = P_R(1,1);
                        end
                        if whichcurrents == 2,
                            P_R_hat = P_R(2,2);
                        end
                        if s_opts.include_ult_radiation_losses
                            P = P + P_R_hat;
                            if s_opts.save_noise_contributions
                                P_R_set(counter:counter+counterincrement) = diag(P_R_hat);
                            end
                        end
                        
                        %                     P_ult(counter_main:counter_main+counterincrement) = diag(P);
                        %                     P_L_modes(counter_main:counter_main+counterincrement) = diag(T_hat*P_L*T_hat');
                        %                     P_A_modes(counter_main:counter_main+counterincrement) = diag(P_A_hat);
                        %                     P_R_modes(counter_main:counter_main+counterincrement) = diag(P_R_hat);
                        %
                        %                     counter_main = counter_main + 1 + counterincrement;
                    end
                    % save P to be used for following voxels % RL 13 Jan 2015
                    phi_matrix_set(:,counter:counter+counterincrement) = P;
                    P_L_matrix_set(:,2*counter_modes-1:2*counter_modes) = P_L;
                    if s_opts.include_ult_conductor_losses || s_opts.include_coil_conductor_losses
                        P_A_matrix_set(2*counter_modes-1:2*counter_modes,1) = diag(P_A_hat);
                    end
                    if s_opts.include_ult_radiation_losses || s_opts.include_coil_radiation_losses
                        P_R_matrix_set(2*counter_modes-1:2*counter_modes,1) = diag(P_R_hat);
                    end
                    
                else % Starting from the second voxel, P is not calculated % RL 13 Jan 2015
                    P = phi_matrix_set(:,counter:counter+counterincrement);
                    P_L = P_L_matrix_set(:,2*counter_modes-1:2*counter_modes);
                    if s_opts.include_ult_conductor_losses || s_opts.include_coil_conductor_losses
                        P_A_hat = diag(P_A_matrix_set(2*counter_modes-1:2*counter_modes,1),0);
                    end
                    if s_opts.include_ult_radiation_losses || s_opts.include_coil_radiation_losses
                        P_R_hat = diag(P_R_matrix_set(2*counter_modes-1:2*counter_modes,1),0);
                    end
                end
                % ------------------------------------------------------- %
                X_P = X.'*inv(P.'); % RL: actually P.' = P because it is a diagonal matrix
                % add X_P and X submatrices to full X_P_full and X_full matrices
                X_P_full(:,counter:counter+counterincrement) = X_P;
                X_full(counter:counter+counterincrement,:) = conj(X); % NB: if whichcurrents = [1 2], then X has the contributions of both the electric and magnetic dipole
                
                if s_opts.compute_coilSNR_flag,
                    % Save quantities which will be needed for computing loop coil SNR
                    if whichcurrents == 1,
                        X_modes(counter_modes,:) = X;
                        %                         P_modes(counter_modes) = P;
                        P_modes(counter_modes) = T_hat*P_L*T_hat'; % RL May 7, 2013
                        if s_opts.include_coil_conductor_losses
%                             P_modes(counter_modes) = P_modes(counter_modes) + P_A_hat;
                            P_modes(counter_modes) = P_modes(counter_modes) + P_A_hat(1,1); % RL 13 Jan 2015
                        end
                        if s_opts.include_coil_radiation_losses
%                             P_modes(counter_modes) = P_modes(counter_modes) + P_R_hat;
                            P_modes(counter_modes) = P_modes(counter_modes) +P_R_hat(1,1); % RL 13 Jan 2015
                        end
                        if s_opts.save_emfields_coil_r4_flag || s_opts.save_emfields_coil_r3_flag || s_opts.save_emfields_coil_r2_flag
                            all_E_phase_modes(counter_modes,:,:) = all_E_phase(counter_modes,:,:);
                            all_B_phase_modes(counter_modes,:,:) = all_B_phase(counter_modes,:,:);
                        end
                    else
                        T_hat = [1 0]*T;
                        X = T_hat*S;
                        X_modes(counter_modes,:) = X;
                        %                         P = T_hat*P_L*T_hat';
                        P_modes(counter_modes)= T_hat*P_L*T_hat';
                        P_L_modes(counter_modes) = T_hat*P_L*T_hat';
                        if s_opts.include_coil_conductor_losses
                            P_modes(counter_modes) = P_modes(counter_modes) + P_A_hat(1,1);
                        end
                        if s_opts.include_coil_radiation_losses % ** TO BE DONE **
                            P_modes(counter_modes) = P_modes(counter_modes) + P_R_hat(1,1);
                            P_R_modes(counter_modes) = P_R_hat(1,1);
                        end
                        if s_opts.save_emfields_coil_r4_flag || s_opts.save_emfields_coil_r3_flag || s_opts.save_emfields_coil_r2_flag
                            if snr_radius == 3, % then all aliasing voxels belong to the same region
                                tempmat = [D_m_41; D_n_41];
                                tempmat = repmat(tempmat,1,length(M_x));
                                % Create EM field components
                                E_x = ele_scaling*(T_hat*(tempmat.*[M_x; N_x]));
                                E_y = ele_scaling*(T_hat*(tempmat.*[M_y; N_y]));
                                E_z = ele_scaling*(T_hat*(tempmat.*[M_z; N_z]));
                                B_x = mag_scaling_r4*(T_hat*(tempmat.*[N_x; M_x]));
                                B_y = mag_scaling_r4*(T_hat*(tempmat.*[N_y; M_y]));
                                B_z = mag_scaling_r4*(T_hat*(tempmat.*[N_z; M_z]));
                                clear tempmat
                            else
                                E_x = zeros(1,length(M_x));
                                E_y = zeros(1,length(M_x));
                                E_z = zeros(1,length(M_x));
                                B_x = zeros(1,length(M_x));
                                B_y = zeros(1,length(M_x));
                                B_z = zeros(1,length(M_x));
                                for irho = 1:length(rset)
                                    if rset(irho) < radius3 % REGION 4
                                        E_x(1,irho) = ele_scaling*(T_hat*([D_m_41; D_n_41].*[M_x(irho); N_x(irho)]));
                                        E_y(1,irho) = ele_scaling*(T_hat*([D_m_41; D_n_41].*[M_y(irho); N_y(irho)]));
                                        E_z(1,irho) = ele_scaling*(T_hat*([D_m_41; D_n_41].*[M_z(irho); N_z(irho)]));
                                        B_x(1,irho) = mag_scaling_r4*(T_hat*([D_m_41; D_n_41].*[N_x(irho); M_x(irho)]));
                                        B_y(1,irho) = mag_scaling_r4*(T_hat*([D_m_41; D_n_41].*[N_y(irho); M_y(irho)]));
                                        B_z(1,irho) = mag_scaling_r4*(T_hat*([D_m_41; D_n_41].*[N_z(irho); M_z(irho)]));
                                    elseif (rset(irho) < radius2) && (rset(irho) >= radius3) % REGION 3
                                        E_x(1,irho) = ele_scaling*(T_hat*[ (B_m_31*M3_x(irho) + D_m_31*M_x(irho)); (B_n_31*N3_x(irho) + D_n_31*N_x(irho)) ]);
                                        E_y(1,irho) = ele_scaling*(T_hat*[ (B_m_31*M3_y(irho) + D_m_31*M_y(irho)); (B_n_31*N3_y(irho) + D_n_31*N_y(irho)) ]);
                                        E_z(1,irho) = ele_scaling*(T_hat*[ (B_m_31*M3_z(irho) + D_m_31*M_z(irho)); (B_n_31*N3_z(irho) + D_n_31*N_z(irho)) ]);
                                        B_x(1,irho) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_x(irho) + D_m_31*N_x(irho)); (B_n_31*M3_x(irho) + D_n_31*M_x(irho)) ]);
                                        B_y(1,irho) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_y(irho) + D_m_31*N_y(irho)); (B_n_31*M3_y(irho) + D_n_31*M_y(irho)) ]);
                                        B_z(1,irho) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_z(irho) + D_m_31*N_z(irho)); (B_n_31*M3_z(irho) + D_n_31*M_z(irho)) ]);
                                    elseif (rset(irho) < radius1) && (rset(irho) >= radius2) % REGION 2
                                        E_x(1,irho) = ele_scaling*(T_hat*[ (B_m_21*M3_x(irho) + D_m_21*M_x(irho)); (B_n_21*N3_x(irho) + D_n_21*N_x(irho)) ]);
                                        E_y(1,irho) = ele_scaling*(T_hat*[ (B_m_21*M3_y(irho) + D_m_21*M_y(irho)); (B_n_21*N3_y(irho) + D_n_21*N_y(irho)) ]);
                                        E_z(1,irho) = ele_scaling*(T_hat*[ (B_m_21*M3_z(irho) + D_m_21*M_z(irho)); (B_n_21*N3_z(irho) + D_n_21*N_z(irho)) ]);
                                        B_x(1,irho) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_x(irho) + D_m_21*N_x(irho)); (B_n_21*M3_x(irho) + D_n_21*M_x(irho)) ]);
                                        B_y(1,irho) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_y(irho) + D_m_21*N_y(irho)); (B_n_21*M3_y(irho) + D_n_21*M_y(irho)) ]);
                                        B_z(1,irho) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_z(irho) + D_m_21*N_z(irho)); (B_n_21*M3_z(irho) + D_n_21*M_z(irho)) ]);
                                    else % REGION 1 (outside the FOV, calculated separately from SNR) 
                                        E_x(1,irho) = 0;
                                        E_y(1,irho) = 0;
                                        E_z(1,irho) = 0;
                                        B_x(1,irho) = 0;
                                        B_y(1,irho) = 0;
                                        B_z(1,irho) = 0;
                                    end
                                end
                            end
                            all_E_phase_modes(counter_modes,:,:) = reshape([E_x E_y E_z],[1 size(xset,2) 3]);
                            all_B_phase_modes(counter_modes,:,:) = reshape([B_x B_y B_z],[1 size(xset,2) 3]);
                        end
                    end
                end
                counter = counter + 1 + counterincrement;
                counter_modes = counter_modes + 1;
            end % end m loop
        end % end l loop
        
        if s_opts.save_emfields_ult_r4_flag || s_opts.save_emfields_ult_r3_flag || s_opts.save_emfields_ult_r2_flag
            % Update the matrix with the modes' EM fields (3 components) at each position
            efield_ult_set(:,ind_1_set,ind_2_set,:) = reshape(all_E_phase,[numbasis length(ind_1_set) length(ind_2_set) 3]);
            bfield_ult_set(:,ind_1_set,ind_2_set,:) = reshape(all_B_phase,[numbasis length(ind_1_set) length(ind_2_set) 3]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ultimate intrinsic SNR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % calculate encoding matrix product
        XPX = X_P_full*X_full;
        %**ADD TEST FOR NANS**
        XPX_inv = inv(XPX);
        % ---------------------------------------------- %
        % calculate SNR denominator for accelerated case
        snr_denom = sqrt(diag(XPX_inv)).';
        snr_ult(ind_1_set,ind_2_set) = reshape(snr_num./snr_denom,[length(ind_1_set) length(ind_2_set)]);
        % ------------------------------------------------ %
        % calculate SNR denominator for unaccelerated case
        snr_denom_unaccel = sqrt(1./diag(XPX)).';
        g_ult(ind_1_set,ind_2_set) = reshape(snr_denom./snr_denom_unaccel,[length(ind_1_set) length(ind_2_set)]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Circular surface coil array SNR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if s_opts.compute_coilSNR_flag
            
            % initialize quantities
            S_coil = zeros(ncoils,length(ind_1_set)*length(ind_2_set));
            psi_coil = zeros(ncoils);
            psi_coil_PR = zeros(ncoils);
            psi_coil_PL = zeros(ncoils);
            
            W_coil = zeros([((lmax + 1)^2 - 1) ncoils]);
            
            
            % define quantities constant for all coils
            costheta_coil_z = coil_offsets(1)/currentradius;
            theta_coil_z = acos(costheta_coil_z);
            %             W_coil_z_norm = (-i*2*pi*coil_radii(1)/(outer_radius^2));
%             W_coil_z_norm = -2*pi*coil_radii(1);
%             
%             % NB: for W_coil_z_norm and costheta_coil_z we assume all coils have the same radius,
%             % otherwise coil_radii(1) and coil_offsets(1) must be replaced with the radius and
%             % the offset of the coil along the z-axis
            
            W_coil_z = zeros(lmax, 1);

            for icoil = 1:ncoils,
                % Compute current weights
                theta_coil_i = rot_coil(icoil,1);
                phi_coil_i = rot_coil(icoil,2);
                costheta_coil_i = cos(theta_coil_i);  %+ coil_offsets(icoil)/outer_radius; 
                
                W_coil_i = zeros([((lmax + 1)^2 - 1) 1]);
                counter_coil_i = 1;
                for l = 1:lmax,
                
                    W_coil_z_norm = 1i*2*pi*sqrt(l/(l+1))/currentradius; % RL 9/9/2013
                                                            
                    if icoil == 1,
                        %**-- calculate the weights for the coil along the z-axis --**
                        legendrefunctions_z = legendre(l,costheta_coil_z,'sch'); % # row is m = 0,...l ; # col is R
                        legendrefunctions_lminus1_z = legendre(l-1,costheta_coil_z,'sch');
                        legendrefunctions_lminus1_z = [legendrefunctions_lminus1_z; zeros(size(costheta_coil_z))];
                        % spherical harmonics with m = 0
                        Y_l_0_z = sqrt((2*l + 1)/(4*pi))*legendrefunctions_z(1,:);
                        Y_lminus1_0_z = sqrt((2*l - 1)/(4*pi))*legendrefunctions_lminus1_z(1,:);

%                         W_coil_z(l) = W_coil_z_norm*(1/(l+1))*(Y_l_0_z*cot(theta_coil_z) - Y_lminus1_0_z*csc(theta_coil_z)*sqrt((2*l + 1)/(2*l - 1)));
                        W_coil_z(l) = W_coil_z_norm*(Y_l_0_z*costheta_coil_z - Y_lminus1_0_z*sqrt((2*l + 1)/(2*l - 1))); % RL 9/9/2013
                        
                        %**----**
                    end
                    rot_coil_norm_i = sqrt(4*pi/(2*l + 1));
                    legendrenorm_i = sqrt((2*l + 1)/(4*pi));
                    legendrefunctions_i = legendre(l,costheta_coil_i,'sch'); % # row is m = 0,...l ; # col is R
                    legendrefunctions_lminus1_i = legendre(l-1,costheta_coil_i,'sch');
                    legendrefunctions_lminus1_i = [legendrefunctions_lminus1_i; zeros(size(costheta_coil_i))];
                    for m = -l:l,
                        if m > 0
                            Y_l_m_i = ((-1)^m)*legendrenorm_i*(1/sqrt(2))*legendrefunctions_i((m+1),:).*exp(1i*m*phi_coil_i);
                        elseif m == 0
                            Y_l_m_i = legendrenorm_i*legendrefunctions_i(1,:);
                        else % using Y_(l,-m) = (-1)^m*conj[Y_(l,m)]
                            Y_l_m_i = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm_i*(1/sqrt(2))*legendrefunctions_i((abs(m)+1),:).*exp(i*abs(m)*phi_coil_i));
                        end

                        W_coil_i(counter_coil_i) = rot_coil_norm_i*conj(Y_l_m_i)*W_coil_z(l);
%                         W_coil_i(counter_coil_i) = ((-1)^m)*rot_coil_norm_i*conj(Y_l_m_i)*W_coil_z(l); % TEST

%                         W_coil(counter_coil_i,icoil) = W_coil_i(counter_coil_i); % RL 16 Oct 2014 (THIS IS INCORRECT, BUT FIXES PROBLEM WITH THE DISPLAY OF COIL CURRENTS)
                        counter_coil_i = counter_coil_i + 1;
                    end % end m loop
                end % end l loop
                
%                 W_coil_i_bis = zeros(size(W_coil_i));
%                 for ll = 1:lmax
%                     startindex = (ll^2 - 1) + 1;
%                     endindex = ((ll+1)^2 - 1);
%                     %                 disp(['l = ' num2str(ll) ', start: ' num2str(startindex) ', end: ' num2str(endindex)]);
%                     for mm = -ll:ll
%                         if mm < 0
%                             index2 = endindex - (ll-abs(mm));
%                             oldidx1 = startindex + abs(abs(mm)-ll);
%                             %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index2)])
%                             W_coil_i_bis(index2) = W_coil_i(oldidx1);
%                         elseif mm > 0
%                             index1 = startindex + abs(mm-ll);
%                             oldidx2 = endindex - (ll-abs(mm));
%                             %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1)])
%                             W_coil_i_bis(index1) = W_coil_i(oldidx2);
%                         elseif mm == 0
%                             index1 = startindex + ll;
%                             %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1)])
%                             W_coil_i_bis(index1) = W_coil_i(index1);
%                         end
%                     end
%                 end
%                 W_coil_i = W_coil_i_bis;
                if s_opts.reverse_m_flag
                    W_coil_i = reverse_m_index(W_coil_i,lmax,1); % RL 25 June 2017: now it is a function
                end
                W_coil(:,icoil) = W_coil_i; % RL 16 Oct 2014 (THIS IS CORRECT, BUT CAUSES A PROBLEM WITH THE DISPLAY OF COIL CURRENTS)
                                
                W_coil_i = W_coil_i(:);
%                 psi_coil(icoil,icoil) = real((abs(W_coil_i.').^2)*P_modes); % check size . should be psi_coil(1,icoil)
                psi_coil(icoil,icoil) = real(W_coil_i.'*(P_modes.*conj(W_coil_i))); % check size . should be psi_coil(1,icoil)
                
%                 psi_coil_PR(1,icoil) = real((abs(W_coil_i.').^2)*P_R_modes); % verifying graph in Keltner 1991
                psi_coil_PR(1,icoil) = real(W_coil_i.'*(P_R_modes.*conj(W_coil_i)));
                %                 psi_coil_PR(1,icoil) = (P_R_modes')*P_R_modes; % verifying graph in Keltner 1991
%                 psi_coil_PL(1,icoil)=real((abs(W_coil_i.').^2)*P_L_modes);
                psi_coil_PL(1,icoil)=real(W_coil_i.'*(P_L_modes.*conj(W_coil_i)));
                PRdivPL=psi_coil_PR./psi_coil_PL;
                S_coil(icoil,:) = W_coil_i.'*X_modes;
%                 all_coil_sens(ind_1,ind_2) = S_coil; % TEMP for revision
                
                if s_opts.save_emfields_coil_r4_flag || s_opts.save_emfields_coil_r3_flag || s_opts.save_emfields_coil_r2_flag
                    
                    for ialiased = 1:length(ind_1_set)*length(ind_2_set)
                        temp_field = (W_coil_i.')*squeeze(all_E_phase_modes(:,ialiased,:));
                        all_E_phase_coil(icoil,ialiased,:) = temp_field;
                        temp_field_2 = (W_coil_i.')*squeeze(all_B_phase_modes(:,ialiased,:));
                        all_B_phase_coil(icoil,ialiased,:) = temp_field_2;
                    end
                    efield_coil_set(icoil,ind_1_set,ind_2_set,:) = reshape(all_E_phase_coil(icoil,:,:),[1 length(ind_1_set) length(ind_2_set) 3]);
                    bfield_coil_set(icoil,ind_1_set,ind_2_set,:) = reshape(all_B_phase_coil(icoil,:,:),[1 length(ind_1_set) length(ind_2_set) 3]); 
                end
                                
                for jcoil = icoil+1:ncoils,

                    theta_coil_j = rot_coil(jcoil,1);
                    phi_coil_j = rot_coil(jcoil,2);
                    costheta_coil_j = cos(theta_coil_j); % coil_offsets(jcoil)/outer_radius;
                    
                    W_coil_j = zeros([((lmax + 1)^2 - 1) 1]);
                    counter_coil_j = 1;
                    for l = 1:lmax
                        rot_coil_norm_j = sqrt(4*pi/(2*l + 1));
                        legendrenorm_j = sqrt((2*l + 1)/(4*pi));
                        legendrefunctions_j = legendre(l,costheta_coil_j,'sch'); % # row is m = 0,...l ; # col is R
                        legendrefunctions_lminus1_j = legendre(l-1,costheta_coil_j,'sch');
                        legendrefunctions_lminus1_j = [legendrefunctions_lminus1_j; zeros(size(costheta_coil_j))];
                        for m = -l:l,
                            %                         for pino = -l:l,
                            %                             m = -pino;
                            if m > 0
                                Y_l_m_j = ((-1)^m)*legendrenorm_j*(1/sqrt(2))*legendrefunctions_j((m+1),:).*exp(i*m*phi_coil_j);
                            elseif m == 0
                                Y_l_m_j = legendrenorm_j*legendrefunctions_j(1,:);
                            else % using Y_(l,-m) = (-1)^m*conj[Y_(l,m)]
                                Y_l_m_j = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm_j*(1/sqrt(2))*legendrefunctions_j((abs(m)+1),:).*exp(i*abs(m)*phi_coil_j));
                            end
                            W_coil_j(counter_coil_j) = rot_coil_norm_j*conj(Y_l_m_j)*W_coil_z(l);
%                             W_coil_j(counter_coil_j) = ((-1)^m)*rot_coil_norm_j*conj(Y_l_m_j)*W_coil_z(l); % TEST
                            counter_coil_j = counter_coil_j + 1;
                        end % end m loop
                    end % end l loop
                    
%                     W_coil_j_bis = zeros(size(W_coil_j));
%                     for ll = 1:lmax
%                         startindex = (ll^2 - 1) + 1;
%                         endindex = ((ll+1)^2 - 1);
%                         %                 disp(['l = ' num2str(ll) ', start: ' num2str(startindex) ', end: ' num2str(endindex)]);
%                         for mm = -ll:ll
%                             if mm < 0
%                                 index2 = endindex - (ll-abs(mm));
%                                 oldidx1 = startindex + abs(abs(mm)-ll);
%                                 %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index2)])
%                                 W_coil_j_bis(index2) = W_coil_j(oldidx1);
%                             elseif mm > 0
%                                 index1 = startindex + abs(mm-ll);
%                                 oldidx2 = endindex - (ll-abs(mm));
%                                 %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1)])
%                                 W_coil_j_bis(index1) = W_coil_j(oldidx2);
%                             elseif mm == 0
%                                 index1 = startindex + ll;
%                                 %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1)])
%                                 W_coil_j_bis(index1) = W_coil_j(index1);
%                             end
%                         end
%                     end
%                     W_coil_j = W_coil_j_bis;
                    if s_opts.reverse_m_flag
                        W_coil_j = reverse_m_index(W_coil_j,lmax,1); % RL 25 June 2017: now it is a function
                    end
                    
                    W_coil_j = W_coil_j(:);
                    psi_coil(icoil,jcoil) = W_coil_i.'*(P_modes.*conj(W_coil_j));
                    psi_coil(jcoil,icoil) = conj(psi_coil(icoil,jcoil));
                end % end jcoil loop
            end % end icoil loop
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate encoding matrix product
            SpsiS_coil = S_coil'*inv(psi_coil)*S_coil;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate SNR denominator for accelerated case
            if s_opts.include_psi_coil
                snr_denom_coil = sqrt(diag(inv(SpsiS_coil))).';
            else
                snr_denom_coil = sqrt(diag( pinv(S_coil)*psi_coil*(pinv(S_coil)')));
                snr_denom_coil = snr_denom_coil.';
            end           
            snr_coil(ind_1_set,ind_2_set) = reshape(snr_num./snr_denom_coil,[length(ind_1_set) length(ind_2_set)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate SNR denominator for unaccelerated case
            % NB: for g-factor the include_psi_coil must be 1
            if s_opts.include_psi_coil
                snr_denom_unaccel_coil = sqrt(1./diag(SpsiS_coil)).';
            else
                snr_denom_unaccel_coil = sqrt(1./diag( pinv(S_coil)*psi_coil*(pinv(S_coil)'))).';
            end
            g_coil(ind_1_set,ind_2_set) = reshape(snr_denom_coil./snr_denom_unaccel_coil,[length(ind_1_set) length(ind_2_set)]);
        else
            snr_coil = NaN;
            g_coil = NaN;
            W_coil = [];
        end
                
        if s_opts.save_ult_current_weights_flag || s_opts.compute_ult_current_pattern_flag || s_opts.compute_ult_field_patterns...
                || s_opts.compute_coil_current_pattern_flag || s_opts.save_coil_current_weights_flag || s_opts.compute_coil_field_patterns,
            % --------------------------------------------------- %
            % compute optimal weights for current basis elements
            wts = (XPX_inv*X_P_full).';
            
            % *** NB: need to change m with -m to calculate the weights *****
            
%             wts_bis = zeros(size(wts));
%             for ll = 1:lmax
% %                 startindex = 2*(ll^2 - 1) + 1;
% %                 endindex = 2*((ll+1)^2 - 1);
%                 % RL July 24, 2015: fixed to work for only one current type
%                 startindex = length(whichcurrents)*(ll^2 - 1) + 1;
%                 endindex = length(whichcurrents)*((ll+1)^2 - 1);
% %                 disp(['l = ' num2str(ll) ', start: ' num2str(startindex) ', end: ' num2str(endindex)]);
%                 for mm = -ll:ll
%                     if mm < 0
%                         index2 = endindex - 2*(ll-abs(mm));
%                         index1 = index2 - 1;
%                         oldidx1 = startindex + 2*abs(abs(mm)-ll);
%                         oldidx2 = 1 + startindex + 2*abs(abs(mm)-ll);
% %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
%                         wts_bis(index1:index2) = wts(oldidx1:oldidx2);
%                     elseif mm > 0
%                         index1 = startindex + 2*abs(mm-ll);
%                         index2 = index1 + 1;
%                         oldidx1 = endindex - 2*(ll-abs(mm)) - 1;
%                         oldidx2 = endindex - 2*(ll-abs(mm));
% %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
%                         wts_bis(index1:index2) = wts(oldidx1:oldidx2);
%                     elseif mm == 0
%                         index1 = startindex + ll*2;
%                         index2 = index1 + 1;
% %                         disp(['    m = ' num2str(mm) ', index = ' num2str(index1) ' and ' num2str(index2)])
%                         wts_bis(index1:index2) = wts(index1:index2);
%                     end
%                 end
%             end
%             wts = wts_bis; % RL Sept 2, 2016 (With this the ideal current patterns are consistent with the 3D plot coordinates)
            
            if s_opts.reverse_m_flag
                wts = reverse_m_index(wts,lmax,whichcurrents,2); % RL 25 June 2017: now it is a function
            end

            wts = reshape(wts,[size(wts,1) length(ind_1_set) length(ind_2_set)]); % size(wts,1) gives the total number of modes

            if s_opts.save_ult_current_weights_flag
                % save weights
                for i_1 = ind_1_set,
                    test_1 = (whichvoxelweights(:,1) == i_1);
                    if any(test_1),
                        for i_2 = ind_2_set,
                            test_2 = (whichvoxelweights(:,2) == i_2);
                            if any(test_1 & test_2),
                                weights(:,:,find(test_1 & test_2)) = ...
                                    permute(...
                                    reshape(...
                                    wts(:,find(ind_1_set == i_1),find(ind_2_set == i_2)),...
                                    [length(whichcurrents) ((lmax + 1)^2 - 1)]),...
                                    [2 1]);
                            end
                        end
                    end
                end
            end
            
            if s_opts.compute_coilSNR_flag && (s_opts.compute_coil_current_pattern_flag || s_opts.save_coil_current_weights_flag || s_opts.compute_coil_field_patterns),
                % compute optimal weights for cylindrical window coils
                wts_coil = (inv(SpsiS_coil)*S_coil'*inv(psi_coil)).';
                if 0,
                    wts_coil = zeros(size(wts_coil));
                    wts_coil(1,:) = 1;
                end
                wts_coil = reshape(wts_coil,[size(wts_coil,1) length(ind_1_set) length(ind_2_set)]);
                if s_opts.save_coil_current_weights_flag
                    for i_1 = ind_1_set,
                        test_1 = (whichvoxelweights(:,1) == i_1);
                        if any(test_1),
                            for i_2 = ind_2_set,
                                test_2 = (whichvoxelweights(:,2) == i_2);
                                if any(test_1 & test_2),
                                    weights_coil(:,find(test_1 & test_2)) = ...
                                        wts_coil(:,find(ind_1_set == i_1),find(ind_2_set == i_2));
                                end
                            end
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if s_opts.compute_ult_current_pattern_flag || s_opts.compute_ult_field_patterns
                % compute current patterns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cptnorm = currentradius;
                %** Note that normalization of phi and theta components of current distribution, and overall normalization, remain to be confirmed.
                for i_1 = ind_1_set,
                    test_1 = (whichvoxels(:,1) == ind_1);
                    if any(test_1),
                        for i_2 = ind_2_set,
                            test_2 = (whichvoxels(:,2) == i_2);
                            if any(test_1 & test_2),
                                disp('    computing current patterns...')
                                if s_opts.compute_ult_field_patterns
                                    disp('    and field patterns...')
                                end
                                if s_opts.compute_patterns_with_stored_weights,
                                    load(weights_file,'currentweights');
                                else
                                    currentweights = ...
                                        permute(...
                                        reshape(...
                                        wts(:,find(ind_1_set == i_1),find(ind_2_set == i_2)),...
                                        [length(whichcurrents) ((lmax + 1)^2 - 1)]),...
                                        [2 1]);
                                end
                                if s_opts.compute_coilSNR_flag && (s_opts.compute_coil_current_pattern_flag || s_opts.compute_coil_field_patterns)
                                    currentweights_coil = wts_coil(:,find(ind_1_set == i_1),find(ind_2_set == i_2));
                                else
                                    currentweights_coil = [];
                                end
                                tic
                                counter_fields = 1;
                                for l = 1:lmax,
                                    disp(['       l = ' num2str(l) '(' num2str(round(toc,1)) ' sec)'])
                                    %                                 for l = 10,
                                    lnorm = sqrt(l*(l+1));
                                    legendrenorm = sqrt((2*l + 1)/(4*pi));
                                    legendrenorm_minus1 = sqrt((2*l - 1)/(4*pi));
                                    if s_opts.compute_ult_current_pattern_flag
                                        legendrefunctions = legendre(l,cos_currenttheta,'sch'); % # row is m = 0,...l ; # col is R
                                        legendrefunctions_lminus1 = zeros([(l + 1) size(cos_currenttheta)]);
                                        legendrefunctions_lminus1(1:l,:,:) = legendre(l-1,cos_currenttheta,'sch');
                                        %                                         legendrefunctions_lminus1 = [legendrefunctions_lminus1; zeros(size(cos_currenttheta))];
                                    end
                                    if s_opts.compute_ult_field_patterns
                                        legendrefunctions_fields = cell(3,1);
                                        legendrefunctions_lminus1_fields = cell(3,1);
                                        for iplane = 1:3
                                            legendrefunctions_fields{iplane} = legendre(l,cosfieldtheta{iplane},'sch'); % # row is m = 0,...l ; # col is R
                                            legendrefunctions_lminus1_fields{iplane} = zeros([(l + 1) size(cosfieldtheta{iplane})]);
                                            legendrefunctions_lminus1_fields{iplane}(1:l,:,:,:) = legendre(l-1,cosfieldtheta{iplane},'sch');
                                        end
                                    end
                                    for m = -l:l,
                                        %%--------------------------%%
                                        %  COMPUTE CURRENT PATTERNS  %
                                        %%--------------------------%%
                                        if s_opts.compute_ult_current_pattern_flag
                                            
                                            lmul = sqrt((2*l + 1)*(l^2 - m^2)/(2*l - 1));
                                            if m > 0
                                                Y_l_m = ((-1)^m)*legendrenorm*(1/sqrt(2))*squeeze(legendrefunctions((m+1),:,:)).*exp(1i*m*currentphi);
                                            elseif m == 0
                                                Y_l_m = legendrenorm*squeeze(legendrefunctions(1,:,:));
                                            else % using Y_(l,-m) = (-1)^m*conj[Y_(l,m)]
                                                Y_l_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm*(1/sqrt(2))*squeeze(legendrefunctions((abs(m)+1),:,:)).*exp(1i*abs(m)*currentphi));
                                            end
                                            if m > 0
                                                Y_lminus1_m = ((-1)^m)*legendrenorm_minus1*(1/sqrt(2))*squeeze(legendrefunctions_lminus1(m+1,:,:)).*exp(1i*m*currentphi);
                                            elseif m == 0
                                                Y_lminus1_m = legendrenorm_minus1*squeeze(legendrefunctions_lminus1(1,:,:));
                                            else
                                                Y_lminus1_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm_minus1*(1/sqrt(2))*squeeze(legendrefunctions_lminus1((abs(m)+1),:,:)).*exp(1i*abs(m)*currentphi));
                                            end
                                            % form and weight current patterns
                                            cpt = zeros([size(currenttheta) 2]);
                                            if s_opts.compute_coil_current_pattern_flag
                                                cpt_coil = cpt;
                                            end
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            % ultimate intrinsic current pattern
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %-- remove NaN at theta=0 --%
                                            [smtheta_1, smtheta_2] = find(abs(currenttheta) < 0.01);
                                            csc_Y_l_m = csc_currenttheta.*Y_l_m;
                                            csc_Y_lminus1_m = csc_currenttheta.*Y_lminus1_m;
                                            if m == -1
                                                csc_Y_l_m(smtheta_1, smtheta_2) = csc_Y_l_m_series_expans(currenttheta(smtheta_1,smtheta_2),currentphi(smtheta_1,smtheta_2),l,m);
                                                csc_Y_lminus1_m(smtheta_1, smtheta_2) = csc_Y_l_m_series_expans(currenttheta(smtheta_1,smtheta_2),currentphi(smtheta_1,smtheta_2),l-1,m);
                                            elseif m == 1
                                                csc_Y_l_m(smtheta_1, smtheta_2) = -conj(csc_Y_l_m_series_expans(currenttheta(smtheta_1,smtheta_2),currentphi(smtheta_1,smtheta_2),l,-m));
                                                csc_Y_lminus1_m(smtheta_1, smtheta_2) = -conj(csc_Y_l_m_series_expans(currenttheta(smtheta_1,smtheta_2),currentphi(smtheta_1,smtheta_2),l-1,-m));
                                            else
                                                csc_Y_l_m(smtheta_1, smtheta_2) = 0;
                                                csc_Y_lminus1_m(smtheta_1, smtheta_2) = 0;
                                            end
                                            %----------------------------%
                                            if whichcurrents == 1,
                                                %                                             cpt(:,:,1) = 1i*m*csc_Y_l_m;
                                                %                                             cpt(:,:,2) = -l*cos_currenttheta.*csc_Y_l_m + lmul*csc_Y_lminus1_m;
                                                cpt(:,:,1) = (1i/lnorm)*1i*m*csc_Y_l_m;
                                                cpt(:,:,2) = (1i/lnorm)*(-l*cos_currenttheta.*csc_Y_l_m + lmul*csc_Y_lminus1_m);
                                                cpt = currentweights( (l^2 -1) + (m + l + 1) )*cpt;
                                            elseif whichcurrents == 2,
                                                %                                             cpt(:,:,1) = l*cos_currenttheta.*csc_Y_l_m - lmul*csc_Y_lminus1_m;
                                                %                                             cpt(:,:,2) = 1i*m*csc_Y_l_m;
                                                cpt(:,:,1) = (1i/lnorm)*(l*cos_currenttheta.*csc_Y_l_m - lmul*csc_Y_lminus1_m);
                                                cpt(:,:,2) = (1i/lnorm)*1i*m*csc_Y_l_m;
                                                cpt = currentweights( (l^2 -1) + (m + l + 1) )*cpt;
                                            else
                                                %                                             cpt(:,:,1) = 1i*m*csc_Y_l_m.*(currentweights((l^2 -1) + (m + l + 1),1)) + ...
                                                %                                                 (l*cos_currenttheta.*csc_Y_l_m - lmul*csc_Y_lminus1_m).*currentweights((l^2 -1) + (m + l + 1),2);
                                                %                                             cpt(:,:,2) = (-l*cos_currenttheta.*csc_Y_l_m + lmul*csc_Y_lminus1_m).*currentweights((l^2 -1) + (m + l + 1),1) + ...
                                                %                                                 1i*m*csc_Y_l_m.*currentweights((l^2 -1) + (m + l + 1),2);
                                                cpt(:,:,1) = (1i/lnorm)*(1i*m*csc_Y_l_m.*(currentweights((l^2 -1) + (m + l + 1),1)) + ...
                                                    (l*cos_currenttheta.*csc_Y_l_m - lmul*csc_Y_lminus1_m).*currentweights((l^2 -1) + (m + l + 1),2));
                                                cpt(:,:,2) = (1i/lnorm)*((-l*cos_currenttheta.*csc_Y_l_m + lmul*csc_Y_lminus1_m).*currentweights((l^2 -1) + (m + l + 1),1) + ...
                                                    1i*m*csc_Y_l_m.*currentweights((l^2 -1) + (m + l + 1),2));
                                                
                                                
                                                %                                             cpt(:,:,1) = i*m*csc_currenttheta.*Y_l_m.*(currentweights((l^2 -1) + (m + l + 1),1)) + ...
                                                %                                                 (l*cot_currenttheta.*Y_l_m - lmul*csc_currenttheta.*Y_lminus1_m).*currentweights((l^2 -1) + (m + l + 1),2)*outer_radius;
                                                %                                             cpt(:,:,2) = (-l*cot_currenttheta.*Y_l_m + lmul*csc_currenttheta.*Y_lminus1_m).*currentweights((l^2 -1) + (m + l + 1),1) + ...
                                                %                                                 i*m*csc_currenttheta.*Y_l_m.*currentweights((l^2 -1) + (m + l + 1),2)*outer_radius;
                                            end
                                            % normalize current densities for discrete display (i.e. multiply by cross-sectional area of normal current pattern voxel faces)
                                            cpt(:,:,1) = cpt(:,:,1)*cptnorm;
                                            cpt(:,:,2) = cpt(:,:,2)*cptnorm;
                                            %  add weighted currents to running pattern buffer
                                            currentpattern(:,:,:,find(test_1 & test_2)) = currentpattern(:,:,:,find(test_1 & test_2)) + cpt;
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            % circular coils current pattern
                                            if s_opts.compute_coil_current_pattern_flag,
                                                cpt_coil(:,:,1) = (1i/lnorm)*(1i*m*csc_Y_l_m);
                                                cpt_coil(:,:,2) = (1i/lnorm)*(-l*cos_currenttheta.*csc_Y_l_m + lmul*csc_Y_lminus1_m);
                                                % normalize current densities for discrete display (i.e. multiply by cross-sectional area of normal current pattern voxel faces)
                                                cpt_coil(:,:,1) = cpt_coil(:,:,1)*cptnorm;
                                                cpt_coil(:,:,2) = cpt_coil(:,:,2)*cptnorm;
                                                %                                             disp(['mode #: ' num2str((l^2 -1) + (m + l + 1))]);
                                                netcurrentweights_coil = W_coil((l^2 -1) + (m + l + 1),:)*currentweights_coil(:);
                                                %                                             netcurrentweights_coil = ((-1)^m)*W_coil((l^2 -1) + (m + l + 1),:)*currentweights_coil(:); %TEST
                                                cpt_coil = netcurrentweights_coil*cpt_coil;
                                                currentpattern_coil(:,:,:,find(test_1 & test_2)) = currentpattern_coil(:,:,:,find(test_1 & test_2)) + cpt_coil;
                                            else
                                                netcurrentweights_coil = [];
                                            end
                                        end % compute_ult_current_pattern_flag
                                        %%-------------------------%%
                                        %  COMPUTE FIELDS PATTERNS  %
                                        %%-------------------------%%
                                        if s_opts.compute_ult_field_patterns
                                                for iplane = 1:3,
                                                    [Efieldpattern{iplane},Bfieldpattern{iplane},Efieldpattern_coil{iplane},Bfieldpattern_coil{iplane}] = ...
                                                        calculate_field_patterns_function(...
                                                        whichcurrents, s_opts,T_matrix_set,l,m,currentweights,W_coil,currentweights_coil,test_1,test_2,counter_fields,...
                                                        radius3,radius2,radius1,k_0_rad_b,mask_r4{iplane},mask_r3{iplane},mask_r2{iplane},mask_r1{iplane},...
                                                        legendrefunctions_fields{iplane}, legendrefunctions_lminus1_fields{iplane},...
                                                        mag_scaling_r1,mag_scaling_r2,mag_scaling_r3,mag_scaling_r4,ele_scaling,...
                                                        k_0_squared,k_0,k_2_squared,k_2,k_3_squared,k_3,k_4_squared,k_4,...
                                                        fieldpatternmatrixsize{iplane},Efieldpattern{iplane},Bfieldpattern{iplane},Efieldpattern_coil{iplane},Bfieldpattern_coil{iplane},...
                                                        fieldrho{iplane},cosfieldtheta{iplane},fieldtheta{iplane},fieldphi{iplane},sinfieldtheta{iplane},sinfieldphi{iplane},...
                                                        cosfieldphi{iplane},cotfieldtheta{iplane},cscfieldtheta{iplane},cos2fieldtheta{iplane});
                                                end
                                        end % compute_ult_field_patterns
                                        counter_fields = counter_fields + 1;
                                    end % m loop
                                end % l loop
                            end
                        end
                    end
                end
            
            end
        end
    end % end loop indices in phase direction
end % end loop indices in frequency direction

clear phi_matrix_set P_L_matrix_set P_A_matrix_set P_R_matrix_set

if 1 % plotting test
myplane = Efieldpattern{1};
myplane_rms = sqrt(abs(myplane(:,:,:,1)).^2 + abs(myplane(:,:,:,2)).^2 + abs(myplane(:,:,:,3)).^2);
figure;
subplot(2,2,1); imshow(squeeze(mask_r4{1}).*myplane_rms,[]);
subplot(2,2,2); imshow(squeeze(mask_r3{1}).*myplane_rms,[]);
subplot(2,2,3); imshow(squeeze(mask_r2{1}).*myplane_rms,[]);
subplot(2,2,4); imshow(squeeze(mask_r2{1}+mask_r3{1}+mask_r4{1}).*myplane_rms,[]);
colormap(parula)
figure; imshow(squeeze(mask_r2{1}+mask_r3{1}+mask_r4{1}).*squeeze(myplane_rms),[]); colormap(parula)
clear myplane myplane_rms

myplane = Efieldpattern{2};
myplane_rms = sqrt(abs(myplane(:,:,:,1)).^2 + abs(myplane(:,:,:,2)).^2 + abs(myplane(:,:,:,3)).^2);
figure;
subplot(2,2,1); imshow(squeeze(mask_r4{2}).*squeeze(myplane_rms),[]);
subplot(2,2,2); imshow(squeeze(mask_r3{2}).*squeeze(myplane_rms),[]);
subplot(2,2,3); imshow(squeeze(mask_r2{2}).*squeeze(myplane_rms),[]);
subplot(2,2,4); imshow(squeeze(mask_r2{2}+mask_r3{2}+mask_r4{2}).*squeeze(myplane_rms),[]);
colormap(parula)
figure; imshow(squeeze(mask_r2{1}+mask_r3{1}+mask_r4{1}).*squeeze(myplane_rms),[]); colormap(parula)
clear myplane myplane_rms

myplane = Efieldpattern{3};
myplane_rms = sqrt(abs(myplane(:,:,:,1)).^2 + abs(myplane(:,:,:,2)).^2 + abs(myplane(:,:,:,3)).^2);
figure;
subplot(2,2,1); imshow(squeeze(mask_r4{3}).*squeeze(myplane_rms),[]);
subplot(2,2,2); imshow(squeeze(mask_r3{3}).*squeeze(myplane_rms),[]);
subplot(2,2,3); imshow(squeeze(mask_r2{3}).*squeeze(myplane_rms),[]);
subplot(2,2,4); imshow(squeeze(mask_r2{3}+mask_r3{3}+mask_r4{3}).*squeeze(myplane_rms),[]);
colormap(parula)
figure; imshow(squeeze(mask_r2{1}+mask_r3{1}+mask_r4{1}).*squeeze(myplane_rms),[]); colormap(parula)
clear myplane myplane_rms
end

% -- EM fields -- %
if s_opts.save_emfields_ult_r4_flag || s_opts.save_emfields_ult_r3_flag || s_opts.save_emfields_ult_r2_flag
        efield_ult_set_fullfov = zeros(numbasis,nf,np,3);
        efield_ult_set_fullfov(:,min(snr_row):max(snr_row),min(snr_col):max(snr_col),:) = efield_ult_set;
        bfield_ult_set_fullfov = zeros(numbasis,nf,np,3);
        bfield_ult_set_fullfov(:,min(snr_row):max(snr_row),min(snr_col):max(snr_col),:) = bfield_ult_set;
end
clear efield_ult_set bfield_ult_set
if s_opts.save_emfields_coil_r4_flag || s_opts.save_emfields_coil_r3_flag || s_opts.save_emfields_coil_r2_flag
        efield_coil_set_fullfov = zeros(ncoils,nf,np,3);
        efield_coil_set_fullfov(:,min(snr_row):max(snr_row),min(snr_col):max(snr_col),:) = efield_coil_set;
        bfield_coil_set_fullfov = zeros(ncoils,nf,np,3);
        bfield_coil_set_fullfov(:,min(snr_row):max(snr_row),min(snr_col):max(snr_col),:) = bfield_coil_set;
end
clear efield_coil_set bfield_coil_set

%%%%-----------------------------------------------------%%%%
% CALCULATE EM FIELD IN REGIONS NOT INCLUDED IN THE SNR FOV %
%%%%-----------------------------------------------------%%%%

% NOTE: if snr_radius = 1, then EM fields for regions 4 (innermost), 3, 2 have already been calculated

[reg2_x,reg2_y] = find(fov_reg2 ~= 0);
% reg2_indx = find(fov_reg2 ~= 0);

if (snr_radius == 2 || snr_radius == 3) && (s_opts.save_emfields_ult_r2_flag || s_opts.save_emfields_coil_r2_flag)
    % calculate EM field for region 2 only, since it has not been calculated yet
    disp('Starting EM field calculation in region 2 ...');
    xset = zeros(1,length(reg2_x));
    yset = zeros(1,length(reg2_x));
    zset = zeros(1,length(reg2_x));
    for ir2 = 1:length(reg2_x)
        xset(1,ir2) = x_fov(reg2_x(ir2),reg2_y(ir2));
        yset(1,ir2) = y_fov(reg2_x(ir2),reg2_y(ir2));
        zset(1,ir2) = z_fov(reg2_x(ir2),reg2_y(ir2));
    end
    
    % Convert in spherical coordinates
    rset = sqrt(xset.^2 + yset.^2 + zset.^2);   % rho
    costhetaset = zset./rset;                   % cos[theta]
    costhetaset(isnan(costhetaset)) = 0;        %** NaN at r=0
    thetaset = acos(costhetaset);               % theta
    phiset = atan2(yset,xset);                  % phi
    
    sinthetaset = sin(thetaset);                % sin[theta]
    sinphiset = sin(phiset);                    % sin[phi]
    cosphiset = cos(phiset);                    % cos[phi]
    
    cotthetaset = cot(thetaset);                % cot[theta] ** Inf at theta=0
    cscthetaset = csc(thetaset);                % csc[theta] ** Inf at theta=0
    cos2thetaset = costhetaset.*costhetaset;    % cos[theta]^2
    
    % unit vectors conversion
    rhat_x = sinthetaset.*cosphiset;
    rhat_y = sinthetaset.*sinphiset;
    rhat_z = costhetaset;
    
    krset = k_2*rset;
    besselnorm = sqrt(pi./(2.*krset));
    
    counter = 1;
    counter_modes = 1;
    
    if s_opts.compute_coilSNR_flag,
        clear all_E_phase_modes all_E_phase_coil all_B_phase_modes all_B_phase_coil
        all_E_phase_modes = zeros(((lmax + 1)^2 - 1),length(reg2_x),3);
        all_B_phase_modes = zeros(((lmax + 1)^2 - 1),length(reg2_x),3);
    end
    
    for l = 1:lmax, % for l=0 the T matrix is empty and there would be no contributions
        %             disp(['      l = ' num2str(l)])
        lnorm = sqrt(l*(l+1));
        legendrenorm = sqrt((2*l + 1)/(4*pi));
        legendrenorm_minus1 = sqrt((2*l - 1)/(4*pi));
        legendrefunctions = legendre(l,costhetaset,'sch'); % # row is m = 0,...l ; # col is R
        if l>0,
            legendrefunctions_lminus1 = legendre(l-1,costhetaset,'sch');
            legendrefunctions_lminus1 = [legendrefunctions_lminus1; zeros(size(costhetaset))];
        end
        
        H_l_1_k_0_rad_b = besselnorm_rad_b*besselh(l+0.5,1,k_0_rad_b);
        H_l_1_k_0_rad_b_prime = -k_0_rad_b*besselnorm_rad_b*besselh(l+1+0.5,1,k_0_rad_b) + (l + 1)*H_l_1_k_0_rad_b;
        
        J_l_kr = besselnorm.*besselj(l+0.5,krset);
        J_l_kr(rset<eps) = (l==0);                            %** only j_0 survives at r=0
        J_lplus1_kr = besselnorm.*besselj(l+1+0.5,krset);
        J_lplus1_kr(rset<eps) = 0;                           %** j_l+1 vanishes at r=0 for l>=0
        J_l_kr_prime = -krset.*J_lplus1_kr + (l + 1)*J_l_kr;
        
        H_l_1_kr = besselnorm.*besselh(l+0.5,1,krset);
        H_l_1_kr_prime = -krset.*besselnorm.*besselh(l+1+0.5,1,krset) + (l + 1)*H_l_1_kr;
        
        [R_H_P_1, R_H_F_1, R_V_P_1, R_V_F_1, T_H_P_1, T_H_F_1, T_V_P_1, T_V_F_1] = ...
            compute_reflection_and_transmission_coef(l,k_0,k_2,g_opts.radius1);
        
        [R_H_P_2, R_H_F_2, R_V_P_2, R_V_F_2, T_H_P_2, T_H_F_2, T_V_P_2, T_V_F_2] = ...
            compute_reflection_and_transmission_coef(l,k_2,k_3,g_opts.radius2);
        
        [R_H_P_3, R_H_F_3, R_V_P_3, R_V_F_3, T_H_P_3, T_H_F_3, T_V_P_3, T_V_F_3] = ...
            compute_reflection_and_transmission_coef(l,k_3,k_4,radius3);
        
        B_m_11 = -( T_H_P_2*(T_H_P_1*R_H_F_1 + T_H_F_1*R_H_F_2) + R_H_F_3*T_H_F_2*(T_H_F_1 + T_H_P_1*R_H_P_2*R_H_F_1) )/...
            ( T_H_P_2*(T_H_P_1 + T_H_F_1*R_H_F_2*R_H_P_1) + R_H_F_3*T_H_F_2*(T_H_P_1*R_H_P_2 + T_H_F_1*R_H_P_1) );
        
        
        B_n_11 = -( T_V_P_2*(T_V_P_1*R_V_F_1 + T_V_F_1*R_V_F_2) + R_V_F_3*T_V_F_2*(T_V_F_1 + T_V_P_1*R_V_P_2*R_V_F_1) )/...
            ( T_V_P_2*(T_V_P_1 + T_V_F_1*R_V_F_2*R_V_P_1) + R_V_F_3*T_V_F_2*(T_V_P_1*R_V_P_2 + T_V_F_1*R_V_P_1) );
        
        B_m_21 = (1/T_H_F_1)*(B_m_11 + R_H_F_1);
        
        B_n_21 = (1/T_V_F_1)*(B_n_11 + R_V_F_1);
        
        D_m_21 = (1/T_H_P_1)*(1+ R_H_P_1*B_m_11);
        
        D_n_21 = (1/T_V_P_1)*(1+ R_V_P_1*B_n_11);
        
        B_m_31 = (1/T_H_F_2)*(B_m_21 + R_H_F_2*D_m_21);
        
        B_n_31 = (1/T_V_F_2)*(B_n_21 + R_V_F_2*D_n_21);
        
        D_m_31 = (1/T_H_P_2)*(R_H_P_2*B_m_21 + D_m_21);
        
        D_n_31 = (1/T_V_P_2)*(R_V_P_2*B_n_21 + D_n_21);
        
        D_m_41 = (1/T_H_P_3)*(R_H_P_3*B_m_31 + D_m_31);
        
        D_n_41 = (1/T_V_P_3)*(R_V_P_3*B_n_31 + D_n_31);
        
        % T matrix does not depend on the expansion parameter m
        Ttemp = [H_l_1_k_0_rad_b                     0
            0            (1/k_0_rad_b)*H_l_1_k_0_rad_b_prime].';
        
        for m = -l:l,
            T = ((-1)^(1-m))*Ttemp;
            
            lmul = sqrt((2*l + 1)*(l^2 - m^2)/(2*l - 1));
            
            if m > 0
                Y_l_m = ((-1)^m)*legendrenorm*(1/sqrt(2))*legendrefunctions((m+1),:).*exp(1i*m*phiset);
            elseif m == 0
                Y_l_m = legendrenorm*legendrefunctions(1,:);
            else % using Y_(l,-m) = (-1)^m*conj[Y_(l,m)]
                Y_l_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm*(1/sqrt(2))*legendrefunctions((abs(m)+1),:).*exp(1i*abs(m)*phiset));
            end
            S_E = zeros(size(Y_l_m));
            S_M = S_E;
            if m > 0
                Y_lminus1_m = ((-1)^m)*legendrenorm_minus1*(1/sqrt(2))*legendrefunctions_lminus1(m+1,:).*exp(1i*m*phiset);
            elseif m == 0
                Y_lminus1_m = legendrenorm_minus1*legendrefunctions_lminus1(1,:);
            else
                Y_lminus1_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm_minus1*(1/sqrt(2))*legendrefunctions_lminus1((abs(m)+1),:).*exp(1i*abs(m)*phiset));
            end
            
            X_x = (1/lnorm).*((-m*cosphiset+1i*l*sinphiset).*cotthetaset.*Y_l_m ...
                - 1i*lmul*cscthetaset.*sinphiset.*Y_lminus1_m);
            X_y = (1/lnorm).*((-m*sinphiset-1i*l*cosphiset).*cotthetaset.*Y_l_m ...
                + 1i*lmul*cscthetaset.*cosphiset.*Y_lminus1_m);
            X_z = (1/lnorm).*m.*Y_l_m;
            
            r_cross_X_x = (1/lnorm).*((m*sinphiset+1i*l*cos2thetaset.*cosphiset).*cscthetaset.*Y_l_m ...
                - 1i*lmul*cotthetaset.*cosphiset.*Y_lminus1_m);
            r_cross_X_y = (1/lnorm).*((-m*cosphiset+1i*l*cos2thetaset.*sinphiset).*cscthetaset.*Y_l_m ...
                - 1i*lmul*cotthetaset.*sinphiset.*Y_lminus1_m);
            r_cross_X_z = (-1i/lnorm).*(l*costhetaset.*Y_l_m ...
                - lmul.*Y_lminus1_m);
            
            M_x = J_l_kr.*X_x;
            M_y = J_l_kr.*X_y;
            M_z = J_l_kr.*X_z;
            
            N_x = (1./krset).*J_l_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_x;
            N_y = (1./krset).*J_l_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_y;
            N_z = (1./krset).*J_l_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_z;
            
            M3_x = H_l_1_kr.*X_x;
            M3_y = H_l_1_kr.*X_y;
            M3_z = H_l_1_kr.*X_z;
            
            N3_x = (1./krset).*H_l_1_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_x;
            N3_y = (1./krset).*H_l_1_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_y;
            N3_z = (1./krset).*H_l_1_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_z;
                        
            T_hat = T;
            if whichcurrents == 1,
                T_hat = [1 0]*T;
            end
            if whichcurrents == 2,
                T_hat = [0 1]*T;
            end
            
            % Create electric field components
            E_x = ele_scaling*(T_hat*[(B_m_21*M3_x + D_m_21*M_x); (B_n_21*N3_x + D_n_21*N_x)]);
            E_y = ele_scaling*(T_hat*[(B_m_21*M3_y + D_m_21*M_y); (B_n_21*N3_y + D_n_21*N_y)]);
            E_z = ele_scaling*(T_hat*[(B_m_21*M3_z + D_m_21*M_z); (B_n_21*N3_z + D_n_21*N_z)]);
            emode_temp = reshape([E_x E_y E_z],[(counterincrement+1) size(xset,2) 3]);
            
            % Create electric field components
            B_x = mag_scaling_r2*(T_hat*[(B_m_21*N3_x + D_m_21*N_x); (B_n_21*M3_x + D_n_21*M_x)]);
            B_y = mag_scaling_r2*(T_hat*[(B_m_21*N3_y + D_m_21*N_y); (B_n_21*M3_y + D_n_21*M_y)]);
            B_z = mag_scaling_r2*(T_hat*[(B_m_21*N3_z + D_m_21*N_z); (B_n_21*M3_z + D_n_21*M_z)]);
            bmode_temp = reshape([B_x B_y B_z],[(counterincrement+1) size(xset,2) 3]);
            
            if s_opts.save_emfields_ult_r2_flag
                for ir2 = 1:length(reg2_x)
                    efield_ult_set_fullfov(counter:counter+counterincrement,reg2_x(ir2),reg2_y(ir2),:) = emode_temp(:,ir2,:);
                    bfield_ult_set_fullfov(counter:counter+counterincrement,reg2_x(ir2),reg2_y(ir2),:) = bmode_temp(:,ir2,:);
                end
            end
            clear emode_temp bmode_temp
            if s_opts.compute_coilSNR_flag && s_opts.save_emfields_coil_r2_flag
                if whichcurrents == 1,
                    all_E_phase_modes(counter_modes,:,:) = all_E_phase(counter_modes,:,:);
                    all_B_phase_modes(counter_modes,:,:) = all_B_phase(counter_modes,:,:);
                else
                    T_hat = [1 0]*T;
                    
                    % Create electric field components
                    E_x = ele_scaling*(T_hat*[(B_m_21*M3_x + D_m_21*M_x); (B_n_21*N3_x + D_n_21*N_x)]);
                    E_y = ele_scaling*(T_hat*[(B_m_21*M3_y + D_m_21*M_y); (B_n_21*N3_y + D_n_21*N_y)]);
                    E_z = ele_scaling*(T_hat*[(B_m_21*M3_z + D_m_21*M_z); (B_n_21*N3_z + D_n_21*N_z)]);
                    all_E_phase_modes(counter_modes,:,:) = reshape([E_x E_y E_z],[1 size(xset,2) 3]);
                    
                    % Create electric field components
                    B_x = mag_scaling_r2*(T_hat*[(B_m_21*N3_x + D_m_21*N_x); (B_n_21*M3_x + D_n_21*M_x)]);
                    B_y = mag_scaling_r2*(T_hat*[(B_m_21*N3_y + D_m_21*N_y); (B_n_21*M3_y + D_n_21*M_y)]);
                    B_z = mag_scaling_r2*(T_hat*[(B_m_21*N3_z + D_m_21*N_z); (B_n_21*M3_z + D_n_21*M_z)]);
                    all_B_phase_modes(counter_modes,:,:) = reshape([B_x B_y B_z],[1 size(xset,2) 3]);
                    
                end
            end
            counter = counter + 1 + counterincrement;
            counter_modes = counter_modes + 1;
        end
    end
    if s_opts.compute_coilSNR_flag && s_opts.save_emfields_coil_r2_flag
        for icoil = 1:ncoils
            W_coil_i = W_coil(:,icoil);
            W_coil_i = W_coil_i(:);
            
            for ir2 = 1:length(reg2_x)
                efield_coil_set_fullfov(icoil,reg2_x(ir2),reg2_y(ir2),:) = (W_coil_i.')*squeeze(all_E_phase_modes(:,ir2,:));
                bfield_coil_set_fullfov(icoil,reg2_x(ir2),reg2_y(ir2),:) = (W_coil_i.')*squeeze(all_B_phase_modes(:,ir2,:));
            end
        end
    end
    disp('*** EM field calculation in region 2 done. ***');
    disp('');
end

[reg3_x,reg3_y] = find(fov_reg3 ~= 0);
% reg3_indx = find(fov_reg3 ~= 0);

if (snr_radius == 3) && (s_opts.save_emfields_ult_r3_flag || s_opts.save_emfields_coil_r3_flag)
    % calculate EM field for region 3 only, since it has not been calculated yet

    disp('Starting EM field calculation in region 3 ...');
    xset = zeros(1,length(reg3_x));
    yset = zeros(1,length(reg3_x));
    zset = zeros(1,length(reg3_x));
    for ir3 = 1:length(reg3_x)
        xset(1,ir3) = x_fov(reg3_x(ir3),reg3_y(ir3));
        yset(1,ir3) = y_fov(reg3_x(ir3),reg3_y(ir3));
        zset(1,ir3) = z_fov(reg3_x(ir3),reg3_y(ir3));
    end
    
    % Convert in spherical coordinates
    rset = sqrt(xset.^2 + yset.^2 + zset.^2);   % rho
    costhetaset = zset./rset;                   % cos[theta]
    costhetaset(isnan(costhetaset)) = 0;        %** NaN at r=0
    thetaset = acos(costhetaset);               % theta
    phiset = atan2(yset,xset);                  % phi
    
    sinthetaset = sin(thetaset);                % sin[theta]
    sinphiset = sin(phiset);                    % sin[phi]
    cosphiset = cos(phiset);                    % cos[phi]
    
    cotthetaset = cot(thetaset);                % cot[theta] ** Inf at theta=0
    cscthetaset = csc(thetaset);                % csc[theta] ** Inf at theta=0
    cos2thetaset = costhetaset.*costhetaset;    % cos[theta]^2
    
    % unit vectors conversion
    rhat_x = sinthetaset.*cosphiset;
    rhat_y = sinthetaset.*sinphiset;
    rhat_z = costhetaset;
    
    krset = k_2*rset;
    besselnorm = sqrt(pi./(2.*krset));
    
    counter = 1;
    counter_modes = 1;
    
    if s_opts.compute_coilSNR_flag,
        clear all_E_phase_modes all_E_phase_coil all_B_phase_modes all_B_phase_coil
        all_E_phase_modes = zeros(((lmax + 1)^2 - 1),length(reg3_x),3);
        all_B_phase_modes = zeros(((lmax + 1)^2 - 1),length(reg3_x),3);
    end
    
    for l = 1:lmax, % for l=0 the T matrix is empty and there would be no contributions
        %             disp(['      l = ' num2str(l)])
        lnorm = sqrt(l*(l+1));
        legendrenorm = sqrt((2*l + 1)/(4*pi));
        legendrenorm_minus1 = sqrt((2*l - 1)/(4*pi));
        legendrefunctions = legendre(l,costhetaset,'sch'); % # row is m = 0,...l ; # col is R
        if l>0,
            legendrefunctions_lminus1 = legendre(l-1,costhetaset,'sch');
            legendrefunctions_lminus1 = [legendrefunctions_lminus1; zeros(size(costhetaset))];
        end
        
        H_l_1_k_0_rad_b = besselnorm_rad_b*besselh(l+0.5,1,k_0_rad_b);
        H_l_1_k_0_rad_b_prime = -k_0_rad_b*besselnorm_rad_b*besselh(l+1+0.5,1,k_0_rad_b) + (l + 1)*H_l_1_k_0_rad_b;
        
        J_l_kr = besselnorm.*besselj(l+0.5,krset);
        J_l_kr(rset<eps) = (l==0);                            %** only j_0 survives at r=0
        J_lplus1_kr = besselnorm.*besselj(l+1+0.5,krset);
        J_lplus1_kr(rset<eps) = 0;                           %** j_l+1 vanishes at r=0 for l>=0
        J_l_kr_prime = -krset.*J_lplus1_kr + (l + 1)*J_l_kr;
        
        H_l_1_kr = besselnorm.*besselh(l+0.5,1,krset);
        H_l_1_kr_prime = -krset.*besselnorm.*besselh(l+1+0.5,1,krset) + (l + 1)*H_l_1_kr;
        
        [R_H_P_1, R_H_F_1, R_V_P_1, R_V_F_1, T_H_P_1, T_H_F_1, T_V_P_1, T_V_F_1] = ...
            compute_reflection_and_transmission_coef(l,k_0,k_2,g_opts.radius1);
        
        [R_H_P_2, R_H_F_2, R_V_P_2, R_V_F_2, T_H_P_2, T_H_F_2, T_V_P_2, T_V_F_2] = ...
            compute_reflection_and_transmission_coef(l,k_2,k_3,g_opts.radius2);
        
        [R_H_P_3, R_H_F_3, R_V_P_3, R_V_F_3, T_H_P_3, T_H_F_3, T_V_P_3, T_V_F_3] = ...
            compute_reflection_and_transmission_coef(l,k_3,k_4,radius3);
        
        B_m_11 = -( T_H_P_2*(T_H_P_1*R_H_F_1 + T_H_F_1*R_H_F_2) + R_H_F_3*T_H_F_2*(T_H_F_1 + T_H_P_1*R_H_P_2*R_H_F_1) )/...
            ( T_H_P_2*(T_H_P_1 + T_H_F_1*R_H_F_2*R_H_P_1) + R_H_F_3*T_H_F_2*(T_H_P_1*R_H_P_2 + T_H_F_1*R_H_P_1) );
        
        
        B_n_11 = -( T_V_P_2*(T_V_P_1*R_V_F_1 + T_V_F_1*R_V_F_2) + R_V_F_3*T_V_F_2*(T_V_F_1 + T_V_P_1*R_V_P_2*R_V_F_1) )/...
            ( T_V_P_2*(T_V_P_1 + T_V_F_1*R_V_F_2*R_V_P_1) + R_V_F_3*T_V_F_2*(T_V_P_1*R_V_P_2 + T_V_F_1*R_V_P_1) );
        
        B_m_21 = (1/T_H_F_1)*(B_m_11 + R_H_F_1);
        
        B_n_21 = (1/T_V_F_1)*(B_n_11 + R_V_F_1);
        
        D_m_21 = (1/T_H_P_1)*(1+ R_H_P_1*B_m_11);
        
        D_n_21 = (1/T_V_P_1)*(1+ R_V_P_1*B_n_11);
        
        B_m_31 = (1/T_H_F_2)*(B_m_21 + R_H_F_2*D_m_21);
        
        B_n_31 = (1/T_V_F_2)*(B_n_21 + R_V_F_2*D_n_21);
        
        D_m_31 = (1/T_H_P_2)*(R_H_P_2*B_m_21 + D_m_21);
        
        D_n_31 = (1/T_V_P_2)*(R_V_P_2*B_n_21 + D_n_21);
        
        D_m_41 = (1/T_H_P_3)*(R_H_P_3*B_m_31 + D_m_31);
        
        D_n_41 = (1/T_V_P_3)*(R_V_P_3*B_n_31 + D_n_31);
        
        % T matrix does not depend on the expansion parameter m
        Ttemp = [H_l_1_k_0_rad_b                     0
            0            (1/k_0_rad_b)*H_l_1_k_0_rad_b_prime].';
        
        for m = -l:l,
            T = ((-1)^(1-m))*Ttemp;
            
            lmul = sqrt((2*l + 1)*(l^2 - m^2)/(2*l - 1));
            
            if m > 0
                Y_l_m = ((-1)^m)*legendrenorm*(1/sqrt(2))*legendrefunctions((m+1),:).*exp(1i*m*phiset);
            elseif m == 0
                Y_l_m = legendrenorm*legendrefunctions(1,:);
            else % using Y_(l,-m) = (-1)^m*conj[Y_(l,m)]
                Y_l_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm*(1/sqrt(2))*legendrefunctions((abs(m)+1),:).*exp(1i*abs(m)*phiset));
            end
            S_E = zeros(size(Y_l_m));
            S_M = S_E;
            if m > 0
                Y_lminus1_m = ((-1)^m)*legendrenorm_minus1*(1/sqrt(2))*legendrefunctions_lminus1(m+1,:).*exp(1i*m*phiset);
            elseif m == 0
                Y_lminus1_m = legendrenorm_minus1*legendrefunctions_lminus1(1,:);
            else
                Y_lminus1_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm_minus1*(1/sqrt(2))*legendrefunctions_lminus1((abs(m)+1),:).*exp(1i*abs(m)*phiset));
            end
            
            X_x = (1/lnorm).*((-m*cosphiset+1i*l*sinphiset).*cotthetaset.*Y_l_m ...
                - 1i*lmul*cscthetaset.*sinphiset.*Y_lminus1_m);
            X_y = (1/lnorm).*((-m*sinphiset-1i*l*cosphiset).*cotthetaset.*Y_l_m ...
                + 1i*lmul*cscthetaset.*cosphiset.*Y_lminus1_m);
            X_z = (1/lnorm).*m.*Y_l_m;
            
            r_cross_X_x = (1/lnorm).*((m*sinphiset+1i*l*cos2thetaset.*cosphiset).*cscthetaset.*Y_l_m ...
                - 1i*lmul*cotthetaset.*cosphiset.*Y_lminus1_m);
            r_cross_X_y = (1/lnorm).*((-m*cosphiset+1i*l*cos2thetaset.*sinphiset).*cscthetaset.*Y_l_m ...
                - 1i*lmul*cotthetaset.*sinphiset.*Y_lminus1_m);
            r_cross_X_z = (-1i/lnorm).*(l*costhetaset.*Y_l_m ...
                - lmul.*Y_lminus1_m);
            
            M_x = J_l_kr.*X_x;
            M_y = J_l_kr.*X_y;
            M_z = J_l_kr.*X_z;
            
            N_x = (1./krset).*J_l_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_x;
            N_y = (1./krset).*J_l_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_y;
            N_z = (1./krset).*J_l_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_z;
            
            M3_x = H_l_1_kr.*X_x;
            M3_y = H_l_1_kr.*X_y;
            M3_z = H_l_1_kr.*X_z;
            
            N3_x = (1./krset).*H_l_1_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_x;
            N3_y = (1./krset).*H_l_1_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_y;
            N3_z = (1./krset).*H_l_1_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_z;
                        
            T_hat = T;
            if whichcurrents == 1,
                T_hat = [1 0]*T;
            end
            if whichcurrents == 2,
                T_hat = [0 1]*T;
            end
            
            % Create electric field components
            E_x = ele_scaling*(T_hat*[(B_m_31*M3_x + D_m_31*M_x); (B_n_31*N3_x + D_n_31*N_x)]);
            E_y = ele_scaling*(T_hat*[(B_m_31*M3_y + D_m_31*M_y); (B_n_31*N3_y + D_n_31*N_y)]);
            E_z = ele_scaling*(T_hat*[(B_m_31*M3_z + D_m_31*M_z); (B_n_31*N3_z + D_n_31*N_z)]);
            emode_temp = reshape([E_x E_y E_z],[(counterincrement+1) size(xset,2) 3]);
            
            % Create electric field components
            B_x = mag_scaling_r3*(T_hat*[(B_m_31*N3_x + D_m_31*N_x); (B_n_31*M3_x + D_n_31*M_x)]);
            B_y = mag_scaling_r3*(T_hat*[(B_m_31*N3_y + D_m_31*N_y); (B_n_31*M3_y + D_n_31*M_y)]);
            B_z = mag_scaling_r3*(T_hat*[(B_m_31*N3_z + D_m_31*N_z); (B_n_31*M3_z + D_n_31*M_z)]);
            bmode_temp = reshape([B_x B_y B_z],[(counterincrement+1) size(xset,2) 3]);
            
            if s_opts.save_emfields_ult_r2_flag
                for ir3 = 1:length(reg3_x)
                    efield_ult_set_fullfov(counter:counter+counterincrement,reg3_x(ir3),reg3_y(ir3),:) = emode_temp(:,ir3,:);
                    bfield_ult_set_fullfov(counter:counter+counterincrement,reg3_x(ir3),reg3_y(ir3),:) = bmode_temp(:,ir3,:);
                end
            end
            clear emode_temp bmode_temp
            if s_opts.compute_coilSNR_flag && s_opts.save_emfields_coil_r2_flag
                if whichcurrents == 1,
                    all_E_phase_modes(counter_modes,:,:) = all_E_phase(counter_modes,:,:);
                    all_B_phase_modes(counter_modes,:,:) = all_B_phase(counter_modes,:,:);
                else
                    T_hat = [1 0]*T;
                    
                    % Create electric field components
                    E_x = ele_scaling*(T_hat*[(B_m_31*M3_x + D_m_31*M_x); (B_n_31*N3_x + D_n_31*N_x)]);
                    E_y = ele_scaling*(T_hat*[(B_m_31*M3_y + D_m_31*M_y); (B_n_31*N3_y + D_n_31*N_y)]);
                    E_z = ele_scaling*(T_hat*[(B_m_31*M3_z + D_m_31*M_z); (B_n_31*N3_z + D_n_31*N_z)]);
                    all_E_phase_modes(counter_modes,:,:) = reshape([E_x E_y E_z],[1 size(xset,2) 3]);
                    
                    % Create electric field components
                    B_x = mag_scaling_r3*(T_hat*[(B_m_31*N3_x + D_m_31*N_x); (B_n_31*M3_x + D_n_31*M_x)]);
                    B_y = mag_scaling_r3*(T_hat*[(B_m_31*N3_y + D_m_31*N_y); (B_n_31*M3_y + D_n_31*M_y)]);
                    B_z = mag_scaling_r3*(T_hat*[(B_m_31*N3_z + D_m_31*N_z); (B_n_31*M3_z + D_n_31*M_z)]);
                    all_B_phase_modes(counter_modes,:,:) = reshape([B_x B_y B_z],[1 size(xset,2) 3]);
                    
                end
            end
            counter = counter + 1 + counterincrement;
            counter_modes = counter_modes + 1;
        end
    end
    if s_opts.compute_coilSNR_flag && s_opts.save_emfields_coil_r2_flag
        for icoil = 1:ncoils
            W_coil_i = W_coil(:,icoil);
            W_coil_i = W_coil_i(:);
            
            for ir3 = 1:length(reg3_x)
                efield_coil_set_fullfov(icoil,reg3_x(ir3),reg3_y(ir3),:) = (W_coil_i.')*squeeze(all_E_phase_modes(:,ir3,:));
                bfield_coil_set_fullfov(icoil,reg3_x(ir3),reg3_y(ir3),:) = (W_coil_i.')*squeeze(all_B_phase_modes(:,ir3,:));
            end
        end
    end    
    
    disp('*** EM field calculation in region 3 done. ***');
    disp('');
end

%%%%-----------------------------------------------%%%%
% CALCULATE EM FIELD IN REGION 1 (OUTSIDE THE OBJECT) %
%%%%-----------------------------------------------%%%%

[reg1_x,reg1_y] = find(fov_reg1 ~= 0);
% reg1_indx = find(fov_reg1 ~= 0);

% @@@@ TO BE DONE @@@@


% if  g_opts.field_pattern % added by MV 9/6/2016 This is for fast field calculation of a complete FOV
%     fast_field_pattern;
% else
%     b1m_max_allplanes_air=0;
%     E_max_allplanes_air=0;
%     efield_ult_set_planes=0;
%     bfield_ult_set_planes=0;
%     fieldx=0;
%     fieldy=0;
%     fieldz=0;
%     mask=0;
%     %     int_E=0;
% end

save(filetosave_ult, 'g_opts','s_opts','nf','np','tempindex_r','tempindex_c','snr_ult', 'g_ult', 'mask_snr','snr_radius','fieldstrength','acceleration','currentradius',...
    'radius1','mask_reg1','fov_reg1',...
    'radius2','tissue_r2','sigma_r2','epsilon_r2','mask_reg2','fov_reg2',...
    'radius3','tissue_r3','sigma_r3','epsilon_r3','mask_reg3','fov_reg3',...
    'tissue_r4','sigma_r4','epsilon_r4','mask_reg4','fov_reg4',...
    'fovf','fovp','x_fov','y_fov','z_fov','lmax','whichcurrents','whichvoxels','whichvoxelweights','sigma_coil','weights','currentpattern','currentphi','currenttheta',...
    'efield_ult_set_fullfov','bfield_ult_set_fullfov','Efieldpattern','Bfieldpattern','P_L_set_r4','P_L_set_r3','P_L_set_r2','P_A_set','P_R_set');

if s_opts.compute_coilSNR_flag
    save(filetosave_coil, 'g_opts','s_opts','nf','np','tempindex_r','tempindex_c','snr_coil', 'g_coil', 'mask_snr','snr_radius','fieldstrength','acceleration','currentradius',...
        'radius1','mask_reg1','fov_reg1',...
        'radius2','tissue_r2','sigma_r2','epsilon_r2','mask_reg2','fov_reg2',...
        'radius3','tissue_r3','sigma_r3','epsilon_r3','mask_reg3','fov_reg3',...
        'tissue_r4','sigma_r4','epsilon_r4','mask_reg4','fov_reg4',...
        'fovf','fovp','x_fov','y_fov','z_fov','lmax','whichcurrents','whichvoxels','whichvoxelweights','sigma_coil','weights_coil','currentpattern_coil','currentphi','currenttheta',...
        'efield_coil_set_fullfov','bfield_coil_set_fullfov','Efieldpattern_coil','Bfieldpattern_coil',...
        'd_coil','ncoils','coil_rotations','coil_radii','coil_offsets','W_coil','psi_coil');
end

output_values = struct(...
                       'g_opts',g_opts,...
                       's_opts',s_opts,...
                       'nf',nf,...
                       'np',np,...
                       'snr_matrix_size',[2*tempindex_r 2*tempindex_c],...
                       'whichcurrents',whichcurrents,...
                       'fieldstrength',fieldstrength,...
                       'acceleration',acceleration,...
                       'snr_ult',snr_ult,...
                       'g_ult',g_ult,...
                       'snr_coil',snr_coil,...
                       'g_coil',g_coil,...
                       'psi_coil',psi_coil,...
                       'mask_snr',mask_snr,...
                       'mask_reg4',mask_reg4,...
                       'mask_reg3',mask_reg3,...
                       'mask_reg2',mask_reg2,...
                       'mask_reg1',mask_reg1,...
                       'fov_reg4',fov_reg4,...
                       'fov_reg3',fov_reg3,...
                       'fov_reg2',fov_reg2,...
                       'fov_reg1',fov_reg1,...
                       'whichvoxels',whichvoxels,...
                       'whichvoxelweights',whichvoxelweights,...
                       'weights',weights,...
                       'weights_coil',weights_coil,...
                       'currentpatternmatrixsize',currentpatternmatrixsize,...
                       'currentpattern',currentpattern,...
                       'currentpattern_coil',currentpattern_coil,...
                       'currentphi',currentphi,...
                       'currenttheta',currenttheta,...
                       'epsilon_r2',epsilon_r2,...
                       'sigma_r2',sigma_r2,...
                       'epsilon_r3',epsilon_r2,...
                       'sigma_r3',sigma_r2,...
                       'epsilon_r4',epsilon_r2,...
                       'sigma_r4',sigma_r2,...
                       'x_fov',x_fov,...
                       'y_fov',y_fov,...
                       'z_fov',z_fov,...
                       'efield_ult_set_fullfov',efield_ult_set_fullfov,...
                       'bfield_ult_set_fullfov',bfield_ult_set_fullfov,...
                       'efield_coil_set_fullfov',efield_coil_set_fullfov,...
                       'bfield_coil_set_fullfov',bfield_coil_set_fullfov,...
                       'Efieldpattern_1',Efieldpattern{1},...
                       'Efieldpattern_2',Efieldpattern{2},...
                       'Efieldpattern_3',Efieldpattern{3},...
                       'Bfieldpattern_1',Bfieldpattern{1},...
                       'Bfieldpattern_2',Bfieldpattern{2},...
                       'Bfieldpattern_3',Bfieldpattern{3},...
                       'Efieldpattern_coil_1',Efieldpattern_coil{1},...
                       'Efieldpattern_coil_2',Efieldpattern_coil{2},...
                       'Efieldpattern_coil_3',Efieldpattern_coil{3},...
                       'Bfieldpattern_coil_1',Bfieldpattern_coil{1},...
                       'Bfieldpattern_coil_2',Bfieldpattern_coil{2},...
                       'Bfieldpattern_coil_3',Bfieldpattern_coil{3},...
                       'PRdivPL',PRdivPL);


disp('****------------------****');
disp('done.');
endtime = toc;
disp(['Total elapsed time = ' num2str(floor(endtime/3600)) 'h ' num2str(floor(mod(endtime,3600)/60)) 'm ' num2str(round(mod(mod(endtime,3600),60))) 's' ]);
disp('****------------------****');


% END dgf_sphere_calc_snr.m

