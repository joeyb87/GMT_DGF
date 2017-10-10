% -------------------------------------------------------------------------
%   DGF_sphere_gmt_batch.m
%
%   Set up Global Maxwell Tomography for a multi-layer spherical phantom
%
%   References:
%
%   Serralles J, Lattanzi R, White J, Sodickson DK, Daniel L and Polymeridis A, 
%   Global Maxwell Tomography: a novel technique for electrical properties mapping 
%   based on MR measurements and volume integral equation formulations;
%   2016 IEEE International Symposium on Antennas and Propagation/
%   USNC-­URSI National Radio Science meeting. 
%   Puerto Rico, 26 June-­1 July 2016, p. TH-A5.2A.4.
%
%
% Riccardo Lattanzi
% August 30, 2017
%
% NOTES:
%       
%       
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%  SET COMMONLY USED SIMULATION LOOP PARAMETERS
% -------------------------------------------------------------------------
%--- frequencies and fields ---%
% whichfrequencies = [1 32 64 96 128 192 256 300 350 400]*1e6;
whichfrequencies = [128]*1e6;
fieldstrength_set = whichfrequencies/42.576e6;
fieldstrength_set = [7]; %[1 3 5 7 9 11];
nfields = length(fieldstrength_set);

matrix_size_set = [65 65 65]; % [x y z]

%--- FOUR-LAYER GEOMETRY AND TISSUE PROPERTIES SETTINGS ---%

mask_radius = 1; % percentage(1 = 100%) of FOV to avoid voxel positions on the circle if needed

%___________________________________%
% ** INNER MOST LAYER - REGION 4 ** %
%-----------------------------------%

% Note: It is possible to specify a set of innermost layer radii
%       (radius3_set) and a set of outermost layer radii (currentradius_set), 
%       but the other two radii must be specify as increments with respect to radius3. 

% AVAILABLE TISSUE TYPES (from http://niremf.ifac.cnr.it/tissprop/)%
% 'Average_Brain' 'Cerebellum' 'Cortical Bone' 'Fat' 'Grey_Matter' 'Skin' 'White_Matter'

radius3 = 0.08; % [m] (innermost radius)
tissue_r4 = []; % use tissue_r4 = [] to manually set tissue properties below

if isempty(tissue_r4),
    % SIEMENS phantom at Bo = 2.8941
%     epsilon_r4 = 80.03;
%     sigma_r4 = 0.084;
    % BRAINO phantom at Bo = 2.8941
%     epsilon_r4 = 60; %81.31;
%     sigma_r4 = 0.45;
    % BASTIEN
    epsilon_r4 = 60; %81.31;
    sigma_r4 = 0.45;
else
    epsilon_r4 = NaN;
    sigma_r4 = NaN;
end

%___________________________________________________%
% **  REGION 3 (radius3 < radius < radius2) ** %
%---------------------------------------------------%

radius2 = radius3 + 0.01; % radius3 < radius < radius2
tissue_r3 = []; % use tissue_r3 = [] to manually set tissue properties below
if isempty(tissue_r3),
    epsilon_r3 = 60;
    sigma_r3 = 0.45;
%     epsilon_r3 = 1;
%     sigma_r3 = 0;
else
    epsilon_r3 = NaN;
    sigma_r3 = NaN;
end

%___________________________________________________%
% **  REGION 2 (radius2 < radius < radius1) ** %
%---------------------------------------------------%

radius1 = radius3 + 0.02; % radius3 < radius < radius2
tissue_r2 = [];
% use tissue_r2 = [] to manually set tissue properties below
if isempty(tissue_r2),
    epsilon_r2 = 60;
    sigma_r2 = 0.45;
%     epsilon_r2 = 1;
%     sigma_r2 = 0;
else
    epsilon_r2 = NaN;
    sigma_r2 = NaN;
end

%________________________________________________%
% ** OUTER MOST LAYER - REGION 1 (always AIR) ** %
%------------------------------------------------%

currentradius = 0.12;     % [m] radius of the surface where the coil/ultimate current density is defined (outermost radius)
nradout = length(currentradius_set);
% **--** NOTE **--** NOTE **--** NOTE **--** NOTE **--** NOTE **--** NOTE **--**
% If outer_radius is > radius1, then it means that there is a layer of air
% (always air) also between radius1 and outer_radius
% **--**       **--**     **--**       **--**     **--**       **--**     **--**

%--- allowed current types ---% (1 = magnetic-dipole, 2 = electric-dipole)
% current_set =      {1   1   [1 2]};
current_set =      [1 2];
% current_set =      {2};
                    
%--- order of the expansion ---%
lmax = 10; % # of modes = 2(lmax + 1)^2

% -------------------------------------------------------------------------
%  SET SIMULATION FLAGS AND OPTIONS
% -------------------------------------------------------------------------

measurement_type = 'b1plus';
% 'b1plus' --> use only B1+ (with relative phases) in the cost function
% 'mrsignal' --> use only the MR signal in the cost function
% 'fullmeasurement' --> use both B1+ (with relative phases) and the MR signal in the cost function

shimming_type = 'none';
% 'none' --> transmit fields are combined with unit weights
% 'cp_mode'--> transmit fields are added constructively at the center
% 'homogeneous' --> transmit fields are combined to maximize homogeneity

compute_current_patterns = 0;  % 1 --> compute current patterns with same combination weights used for the transmit fields

preset_coil_geometry_flag = 0;           % 1 --> coil geometry is loaded from file

%--- saved variables ---%
save_currents_flag = 1;                 % 1 --> save current patterns
save_ult_field_patterns_flag = 0;       % 1 --> save ultimate E and B field patterns (net field, not for each mode, always all layers)
save_coil_field_patterns_flag = 0;      % 1 --> save coil E and B field patterns (net field, not for each coil, always all layers)

%--- plotting ---%
plot_stored_results = 0;                % 1 --> plot results loaded from MAT file (using plot options specified here and not in the stored result file; WARNING: no computation performed if this flag = 1)
plot_gmt_all = 1;                        % 1 --> plot ultimate SNR results
plot_ep = 0;                       % 1 --> plot SNR results for circular surface coils
plot_error = 1;          % 1 --> plot optimal current patterns
plot_currents = 0;           % 1 --> plot modes' combination weights for ultimate SNR
plot_currents_3D = 1;                    % 1 --> display currents pattern on tridimensional surfaces (only if current patterns are plotted) 
plot_ult_field_patterns = 0;            % ##### 1 --> plot ultimate E and B field patterns
plot_coil_field_patterns = 0;           % ##### 1 --> plot coil E and B field patterns
plot_absolutefield_flag = 0;            % ##### 1 --> plot absolute magnitude of fields (corresponding to a time average) in addition to real snapshot
plot_coil_number_flag = 0;              % 1 --> add coil number next to coil geometry
plot_coil_setup_flag = 0;               % 1 --> create figure with sphere, coil geometry and FOV used in the simulation

%--- general plotting options ---%
plot_font_size = 18;                    % ##### FontSize used for some of the plots
allow_windowing_flag = 0;               % ##### 1 --> shows sliding bars to adjust windowing on figure windows (i.e. call to imgui, which can be undone for the active figure by entering imgui in the command window)
save_plots = 0;                         % ##### 1 --> save all figures generated in separate files using export_fig
save_plots_format = 'png';              % ##### available formats: ai, bmp, emf, eps, fig, jpg, pdf, png, tiff
save_plots_resolution = 150;            % ##### specify the resolution of the exported figure

%--- current patterns plotting options ---%
real_part_flag = 1;                     % 1 --> plot the real part of the complex-valued current patterns (0 --> plot imaginary part)
show_axes = 1;                          % 1 --> show Cartesian axes in the 3D plot
axes_color = [0 0 0];                   % 1 --> color of Cartesian axes in the 3D plot
plot_bkgcolor = [0.5 0.5 0.5];          % background color for 3D plots [1 1 1] = white, [0.5 0.5 0.5] = gray, [0 0 0] = black
currentpatternmatrixsize = [52 52];     % number of points at which the current vectors are calculated
whichvoxels = [1 1];                    % coordinates of the voxel for which current patterns are calculated (if [] then all voxels of the FOV)
imagebackgroundflag = 1;                % ##### 1 --> plot current patterns with a background image showing current magnitude
patternphase = 90;                      % ##### phase of current or field pattern plotted (degrees)
currentarrowswidth = 1.1;               % "LineWidth" for quiver plot of current patterns
currentarrowscolor = 'k';               % "Color" for quiver plot of current patterns
spherefacecolor = [0.89 0.89 0.89];     % "FaceColor" color of the sphere in 3D current plots. [0.89 0.89 0.89] = Black&White figures, [1 0.46 0] = Ultimate basis set, [0.9 1 1] = Light blue sphere
sphereedgecolor = 'none';               % "EdgeColor" color of the patches' borders in 3D current plots
spherefacealpha = 1;                    % "FaceAlpha" level of transparency of the patches in 3D current plots
sphereview = [-37.5-180 30];            % azimuth and elevation of the viewpoint in 3D current plots

%--- field patterns options ---%
fieldpatternmatrixsize = [65 65 65];    % ##### number of points at which the E and B field vectors are calculated
% fieldpatternfov = [0.4 0.4 0.4];        % FOV (in meters) for field computations.  If empty, FOV will default to [cylinder diameter cylinder diameter z_fov_scaling*cylinder_length] 
fieldpatternthreeplaneoffset = [0 0 0]; % ##### offset from origin of each of the three orthogonal planes in x, y, or z
fieldpatternshowwhichplanes = [1 2 3];  % ##### on which planes to plot field
fieldpatternlogplots = 0;               % ##### 1 --> plot field patterns on log scale
fieldpatternvolumedisplaystyle = 1;     % ##### 1 --> three-plane style, minus shell [default]; 2 --> montage
simplelabelsflag = 0;                   % ##### 1 --> use simple labels for presentation rather than tracking
reverse_m_flag = 0;                     % 1 --> switch modes order m with -m (Note: the fieldpatterns won't be correct, need to fix it)

%--- EM fields plotting options ---%
scale_factor = 2;                       % new matrix size will be scale_factor*matrix_size
interp_method = 'bilinear';             % interpolation method for scaling the original image (nearest, bilinear, bicubic)
display_range_ult = [0 1.4e-5];         % display range of imshow (for coil fields), specified as a two-element vector [LOW HIGH]; [] --> default
display_range_coil = [0 1.4e-5];        % display range of imshow (for coil fields), specified as a two-element vector [LOW HIGH]; [] --> default
plot_coil_b1plus_phase = 1;             % plot the phase of B1+ for all coils (in addition to magnitude)
plot_coil_b1minus_phase = 1;            % plot the phase of B1- for all coils (in addition to magnitude)
plot_profiles = 1;                      % plot vertical and horizontal profiles
logplot = 1;                            % plot absolute value of fields in log10 scale
                                        % (on Windows are available: 'Ideo3','Ideo5','Cinepak','MSVC','RLE','None')
%--- general ---%
user_label = 'RL';                      % user signature that will be added on some of the plots
warningsoffflag = 1;                    % 1 --> set off some MATLAB warnings

% -------------------------------------------------------------------------
% DEFINE DIRECTORIES PATHS
% -------------------------------------------------------------------------

home_dir     = './';
commondir   = [home_dir 'common_routines/'];
plotdir     = [home_dir 'plot_routines/'];
basissetdir = [home_dir 'results_ultimate/'];
circcoildir = [home_dir 'results_coil/'];
moviedir    = [home_dir 'movies/'];
logdir      = [home_dir 'logfiles/'];
geometrydir = [home_dir 'coil_packing/'];
fovdir      = [home_dir 'fovfiles/'];
tissuedir   = [home_dir 'tissue_properties/'];

addpath(plotdir);
addpath(commondir);
if save_log_file
    logfilename = [logdir 'sphere_log_user_' user_label '_' datestr(now,30) '.txt'];
    diary(logfilename);
end

%--- stored data file names ---%
fov_file = [fovdir 'origin_fov.mat'];
% coil_file = [];
coil_file = [geometrydir '8_coils_mv.mat'];
weights_file = [];  % the format of weights_file shold be [ ((lmax + 1)^2 - 1) length(whichcurrents)]

if preset_fov_geometry_flag % file must contain x_fov,y_fov,z_fov
    plot_regions_fov = 0;
    if fov_file
        load(fov_file);
    else
        [fov_filename,fov_filepath] = uigetfile('*.mat','Select the MAT-file with FOV data');
        fov_file = fullfile(fov_filepath,fov_filename);
        load(fov_file);
    end
    matrix_size = size(x_fov);
    matrix_size_set = repmat(matrix_size,size(acceleration_set,1),1);
end

if plot_stored_results % no computation
    disp('**WARNING: plotting stored results, no computation performed**')
    compute_ultSNR_flag = 0;
    compute_coilSNR_flag = 0;
    compute_ult_current_pattern = 0;
    compute_coil_current_pattern = 0;
    compute_ult_field_patterns = 0;
    compute_coil_field_patterns = 0;
end

% -------------------------------------------------------------------------
%  SET GEOMETRY PARAMETERS
% -------------------------------------------------------------------------

%--- FOV ---%
patientposition = 'headfirst';
patientorientation = 'supine';
sliceorientation = 'transverse';
phasedir = 'LR';
% sliceorientation = 'sagittal';
% phasedir = 'FH';

apoff = 0; lroff = 0; fhoff = 0.0015;
apang = 0; lrang = 0; fhang = 0;
image_plane_offset = [apoff lroff fhoff];
image_plane_orientation = [apang lrang fhang];

%--- circular surface coil geometry ---%
if compute_coilSNR_flag
    if preset_coil_geometry_flag,
        if coil_file
            load(coil_file)
        else
            [coil_file,coil_path] = uigetfile('*.mat','Select the MAT-file with coil array geometry',geometrydir);
            load(fullfile(coil_path,coil_file));
        end
        ncoils = size(coil_rotations,1);
        coil_rotations = 180*coil_rotations/pi; % conversion from rad to deg
                
    else % user defined loop coil data
        coil_rotations = ... % [theta phi]
            [90 0
            90 45
            90 90
            90 135
            90 180
            90 225
            90 270
            90 315];
        coil_rotations = ... % [theta phi]
            [90 0];
        ncoils = size(coil_rotations,1);
%         coil_radius_mul = 0.48;
%                 coil_radius_mul = 0.36;
        coil_radius_mul = 0.36;
    end
    %--- check if all varaiables are there ---%
    if any(isnan(ncoils)) || isempty(coil_rotations) || any(isnan(coil_radius_mul))
        disp('-------------------------------------------------------');
        disp('**ERROR** unspecified circular surface coils variable');
        disp('-------------------------------------------------------');
        pause
    end
else
    ncoils = NaN;
    coil_rotations = [];
end

% -------------------------------------------------------------------------
%  SET EXPERIMENTAL PARAMETERS
% -------------------------------------------------------------------------
sigma_copper = 58e6;             % copper conductivity = 58*10^6 S/m
if include_coil_conductor_losses,
    sigma_coil = sigma_copper;   % [ohm^-1][m^-1] conductivity of the shield
    d_coil = NaN;                % d_coil = NaN --> use skin depth of coil conductor
else
    sigma_coil = Inf;
    d_coil = 0.001;
end

%--- Experimental scaling factors ---%
Vvoxel = 0.002*0.002*0.003;                
NEX = 1;
flipangle = pi/9;
BW = 25.6e3;           % total receive bandwidth
noisefactor = 1.2023;  % noisefigure = 10*log10(noisefactor) = 0.8 [dB]

% -------------------------------------------------------------------------
%  GROUPING SIMULATION SETTINGS AND FLAGS 
% -------------------------------------------------------------------------

actual_size_set = matrix_size_set./acceleration_set;
if ~isempty(whichvoxels) && ((any(whichvoxels(:,1) > actual_size_set(:,1))) || (any(whichvoxels(:,2) > actual_size_set(:,2))))
    whichvoxels = [1 1];
    disp('WARNING: whichvoxel exceeded matrix_size_set./acceleration_set so it has been set to [1 1]');
end

if compute_coilSNR_flag && ~(compute_ultSNR_flag)
    compute_ultSNR_flag = 1;
end
if save_ult_currents_flag && ~(compute_ult_current_pattern_flag)
    compute_ult_current_pattern_flag = 1;
end
if save_coil_currents_flag && ~(compute_coil_current_pattern_flag)
    compute_coil_current_pattern_flag = 1;
end
if compute_coil_current_pattern_flag && ~(compute_ult_current_pattern_flag)
    compute_ult_current_pattern_flag = 1;
end
if save_coil_current_weights_flag && ~(save_ult_current_weights_flag)
    save_ult_current_weights_flag = 1;
end
if compute_patterns_with_stored_weights && ~(compute_ult_field_patterns || compute_coil_field_patterns)
    if ~plot_stored_results
        compute_patterns_with_stored_weights = 0;
        disp('WARNING: compute_patterns_with_stored_weights set to zero because field patterns are not computed');
    end
end
if compute_patterns_with_stored_weights
    % the format of weights_file shold be:
    % [ ((lmax + 1)^2 - 1) length(whichcurrents)]
    if isempty(weights_file)
        [wts_file,wts_path] = uigetfile('*.mat','Select the MAT-file with modes'' ultimate combination weights',geometrydir);
        weights_file = fullfile(wts_path,wts_file);
    end
end

simulation_options = struct(...
    'snr_radius',snr_radius,...
    'lmax',lmax,...
    'include_ult_conductor_losses',include_ult_conductor_losses,...
    'include_coil_conductor_losses',include_coil_conductor_losses,...
    'include_ult_radiation_losses',include_ult_radiation_losses,...
    'include_coil_radiation_losses',include_coil_radiation_losses,...
    'include_circuit_losses',include_circuit_losses,...
    'include_psi_coil',include_psi_coil,...
    'include_experimental_snr_scaling',include_experimental_snr_scaling,...
    'noisefactor',noisefactor,...
    'compute_ultSNR_flag',compute_ultSNR_flag,...
    'compute_coilSNR_flag',compute_coilSNR_flag,...
    'compute_ult_current_pattern_flag',compute_ult_current_pattern_flag,...
    'compute_coil_current_pattern_flag',compute_coil_current_pattern_flag,...
    'compute_ult_field_patterns',compute_ult_field_patterns,...
    'compute_coil_field_patterns',compute_coil_field_patterns,...
    'compute_patterns_with_stored_weights',compute_patterns_with_stored_weights,...
    'weights_file',weights_file,...
    'save_emfields_ult_r4_flag',save_emfields_ult_r4_flag,...
    'save_emfields_coil_r4_flag',save_emfields_coil_r4_flag,...
    'save_emfields_ult_r3_flag',save_emfields_ult_r3_flag,...
    'save_emfields_coil_r3_flag',save_emfields_coil_r3_flag,...
    'save_emfields_ult_r2_flag',save_emfields_ult_r2_flag,...
    'save_emfields_coil_r2_flag',save_emfields_coil_r2_flag,...
    'save_emfields_ult_r1_flag',save_emfields_ult_r1_flag,...
    'save_emfields_coil_r1_flag',save_emfields_coil_r1_flag,...
    'save_ult_current_weights_flag',save_ult_current_weights_flag,...
    'save_coil_current_weights_flag',save_coil_current_weights_flag,...
    'save_ult_currents_flag',save_ult_currents_flag,...
    'save_coil_currents_flag',save_coil_currents_flag,...
    'save_noise_contributions',save_noise_contributions,...
    'save_ult_field_patterns_flag',save_ult_field_patterns_flag,...
    'save_coil_field_patterns_flag',save_coil_field_patterns_flag,...
    'reverse_m_flag',reverse_m_flag,...
    'sigma_coil',sigma_coil,...
    'd_coil',d_coil,...
    'tissue_r4',tissue_r4,...
    'epsilon_r4',epsilon_r4,...
    'sigma_r4',sigma_r4,...
    'tissue_r3',tissue_r3,...
    'epsilon_r3',epsilon_r3,...
    'sigma_r3',sigma_r3,...
    'tissue_r2',tissue_r2,...
    'epsilon_r2',epsilon_r2,...
    'sigma_r2',sigma_r2,...
    'virtual_coils_flag', virtual_coils_flag,...
    'vcoils', vcoils,...
    'poynting_rad_scale',poynting_rad_scale,...
    'npatches',npatches);

path_options = struct(...
    'home_dir',home_dir,...
    'commondir',commondir,...
    'plotdir',plotdir,...
    'basissetdir',basissetdir,...
    'circcoildir',circcoildir,...
    'moviedir',moviedir,...
    'logdir',logdir,...
    'fovdir',fovdir,...
    'tissuedir',tissuedir,...
    'geometrydir',geometrydir);

whichvoxelweights = whichvoxels; % consider making it independent from whichvoxels in the future

if ~plot_stored_results % otherwise these options are needed
    if plot_coilSNR && ~(compute_coilSNR_flag)
        plot_coilSNR = 0;
    end
    if plot_coil_current_patterns && ~(compute_coil_current_pattern)
        plot_coil_current_patterns = 0;
    end
    if plot_coil_current_weights && ~(compute_coilSNR_flag)
        plot_coil_current_weights = 0;
    end
end

plotting_options = struct(...
    'user_label',user_label,...
    'plot_font_size',plot_font_size,...
    'plot_ultSNR',plot_ultSNR,...
    'plot_coilSNR',plot_coilSNR,...
    'plot_ult_current_patterns',plot_ult_current_patterns,...
    'plot_ult_current_weights',plot_ult_current_weights',...
    'plot_coil_current_patterns',plot_coil_current_patterns,...
    'plot_coil_current_weights',plot_coil_current_weights,...
    'plot_current_3D',plot_current_3D,...
    'plot_efields_ult',plot_efields_ult,...
    'plot_efields_coil',plot_efields_coil,...
    'plot_bfields_ult',plot_bfields_ult,...
    'plot_bfields_coil',plot_bfields_coil,...
    'plot_coil_number_flag',plot_coil_number_flag,...
    'plot_coil_setup_flag',plot_coil_setup_flag,...
    'plot_regions_fov',plot_regions_fov,...
    'plot_ult_field_patterns',plot_ult_field_patterns,...
    'plot_coil_field_patterns',plot_coil_field_patterns,...
    'plot_absolutefield_flag',plot_absolutefield_flag,...
    'allow_windowing_flag',allow_windowing_flag,...
    'save_plots',save_plots,...
    'save_plots_format',save_plots_format,...
    'save_plots_resolution',save_plots_resolution,...
    'show_axes',show_axes,...
    'axes_color',axes_color,...
    'movie_ult_current',movie_ult_current,...
    'movie_coil_current',movie_coil_current,...
    'movie_window_scaling',movie_window_scaling,...
    'movie_frames',movie_frames,...
    'movie_fps',movie_fps,...
    'movie_compression',movie_compression,...
    'real_part_flag',real_part_flag,...
    'imagebackgroundflag',imagebackgroundflag,...
    'patternphase',patternphase,...
    'plot_bkgcolor',plot_bkgcolor,...
    'currentpatternmatrixsize',currentpatternmatrixsize,...
    'whichvoxels',whichvoxels,...
    'whichvoxelweights',whichvoxelweights,...
    'currentarrowswidth',currentarrowswidth,...
    'currentarrowscolor',currentarrowscolor,...
    'spherefacecolor',spherefacecolor,...
    'sphereedgecolor',sphereedgecolor,...
    'spherefacealpha',spherefacealpha,...
    'sphereview',sphereview,...
    'scale_factor',scale_factor,...
    'interp_method',interp_method,...
    'display_range_ult',display_range_ult,...
    'display_range_coil', display_range_coil,...
    'plot_coil_b1plus_phase',plot_coil_b1plus_phase,...
    'plot_coil_b1minus_phase',plot_coil_b1minus_phase,...
    'fieldpatternmatrixsize',fieldpatternmatrixsize,...
    'fieldpatternthreeplaneoffset',fieldpatternthreeplaneoffset,...
    'fieldpatternshowwhichplanes',fieldpatternshowwhichplanes,...
    'fieldpatternlogplots',fieldpatternlogplots,...
    'fieldpatternvolumedisplaystyle',fieldpatternvolumedisplaystyle,...
    'simplelabelsflag',simplelabelsflag,...
    'plot_profiles',plot_profiles,...
    'logplot',logplot);

% -------------------------------------------------------------------------
%  RUNNING SNR SIMULATION
% -------------------------------------------------------------------------

if warningsoffflag,
    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:divideByZero')
end

% Error handling
dbstop if error

% -- Ultimate SNR and loop coil SNR--
if compute_ultSNR_flag || compute_coilSNR_flag,
    if ~compute_ultSNR_flag
        simulation_options.compute_ultSNR_flag = 1;
    end
    if compute_coil_current_pattern_flag && ~(compute_coilSNR_flag)
    	disp('** Coil currents can''t be computed as compute_coilSNR_flag is zero')
        simulation_options.compute_coil_current_pattern_flag = 0;
    end
    if save_coil_current_weights_flag && ~(compute_coilSNR_flag)
    	disp('** Coil current weights can''t be computed as compute_coilSNR_flag is zero')
        simulation_options.save_coil_current_weights_flag = 0;
    end
    snr_ult = NaN*zeros([nfields nradout nradin nacc matrix_size_set(nacc,:)]);
    g_ult = snr_ult;
    snr_coil = snr_ult;
    g_coil = g_ult;
    disp('Starting ultimate SNR calculation...');
    for ifield = 1:nfields
        fieldstrength = fieldstrength_set(ifield);
        for iacc = 1:size(acceleration_set,1)   
            acceleration = acceleration_set(iacc,:);
            matrix_size = matrix_size_set(iacc,:);
            %--- experimental SNR scaling ---%
            if include_experimental_snr_scaling
                nacq = prod(matrix_size)/prod(acceleration);
                expsnr_num = Vvoxel*sqrt(nacq*NEX)*sin(flipangle)/(noisefactor*sqrt(BW));
            else
                expsnr_num = 1;
            end
            for iradin = 1:length(radius3_set)
                radius3 = radius3_set(iradin);
                
                radius1 = radius1_set(iradin);
                radius2 = radius2_set(iradin);
                if radius2 > radius1
                    disp('***********************************')
                    disp('    ERROR: Radius 2 > Radius 1     ')
                    disp('***********************************')
                    keyboard
                end
                
                for iradout = 1:length(currentradius_set)    
                    currentradius = currentradius_set(iradout);
                    if radius1 > currentradius
                        radius1 = currentradius;
                    disp('***************************************************')
                    disp('    WARNING: Radius 1 reduced to currentradius     ')
                    disp('***************************************************')
                    end
                    whichcurrents = current_set{iradout};
                    outerfov_radius = outerfov_radius_set(iradout);
                    
%                     fov = [mask_radius*2*outerfov_radius mask_radius*2*outerfov_radius];
%                     fovf = fov(1);
%                     fovp = fov(2);
                    
                    if compute_coilSNR_flag
                        coil_surface_radii = ones(ncoils,1)*currentradius;
                        coil_radii = coil_surface_radii*coil_radius_mul;
                        
%                         coil_radii(1) = coil_surface_radii(1)*coil_radius_mul/2; %TEST TEST TEST
                        
                        coil_offsets = sqrt(coil_surface_radii.^2 - coil_radii.^2);
                    else
                        coil_radii = [];
                        coil_offsets = [];
                    end
                    if compute_coilSNR_flag && plot_circular_coils_flag && (~plot_coil_setup_flag)
                        figure;
                        plot_geometry_sphere(radius3,coil_radii,coil_offsets,coil_rotations,[0.1],[0],[0.1],3,'-',[1 0.46 0],0,plot_coil_number_flag)
                    end
                    
                    geometry_options = struct(...
                            'preset_fov_geometry_flag',preset_fov_geometry_flag,...
                            'fov_file',fov_file,...
                            'matrix_size',matrix_size,...
                            'radius3',radius3,...
                            'radius2',radius2,...
                            'radius1',radius1,...
                            'currentradius',currentradius,...
                            'outerfov_radius',outerfov_radius,...
                            'mask_radius',mask_radius,...
                            'phasedir',phasedir,...
                            'patientposition',patientposition,...
                            'patientorientation',patientorientation,...
                            'sliceorientation',sliceorientation,...
                            'image_plane_offset',image_plane_offset,...
                            'image_plane_orientation',image_plane_orientation,...
                            'ncoils',ncoils,...
                            'coil_rotations',coil_rotations,...
                            'coil_radii',coil_radii,...
                            'coil_offsets',coil_offsets);
                        
                    %------------------------------------------------------
                    [output_values] = dgf_sphere_calc_snr(... 
                        whichcurrents,fieldstrength,acceleration,expsnr_num,simulation_options,geometry_options,...
                        path_options,currentpatternmatrixsize,fieldpatternmatrixsize,fieldpatternthreeplaneoffset,...
                        whichvoxels,whichvoxelweights,plot_poynting_sphere_flag,plot_regions_fov);
                    %------------------------------------------------------
                                        
                    dgf_sphere_plot_snr(output_values,path_options,plotting_options,snr_radius)
                    
                end % iradout loop
            end % iradin loop
        end % iacc loop
    end % ifield loop
    disp('Ultimate SNR calculation done.');
end

%--- PLOTS FROM STORED RESULTS ---%
if plot_stored_results
    if plot_ultSNR || plot_ult_current_patterns || plot_ult_current_weights || plot_efields_ult || plot_bfields_ult || movie_ult_current
        
        [ultSNRfile,ultSNRpath] = uigetfile([ basissetdir 'ultimate_snr_results/*.mat'],'Select the MAT-file with ultSNR results');
        load(fullfile(ultSNRpath,ultSNRfile));
        
        ult_output_values = struct(...
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
            'snr_coil',[],...
            'g_coil',[],...
            'psi_coil',[],...
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
            'weights_coil',[],...
            'weights_file',weights_file,...
            'currentpatternmatrixsize',currentpatternmatrixsize,...
            'currentpattern',currentpattern,...
            'currentpattern_coil',[],...
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
            'efield_coil_set_fullfov',[],...
            'bfield_coil_set_fullfov',[],...
            'PRdivPL',NaN);
        
        if plot_regions_fov
            mu = 4*pi*1e-7;         % permeability of free space [Wb][A^-1][m^-1]
            c = 3e8;                % speed of light [m][s]
            epsilon_0 = 1/(mu*c^2); % permittivity [C][V^-1] = [s^4][A^2][m^-2][kg^-2]
            
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
            title(['Voxel Resolution = ' sprintf('%0.2f',1000*fovf/nf) ' [mm]'],'FontSize',16)
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
        
        plotting_options.plot_coilSNR = 0;
        plotting_options.plot_coil_current_patterns = 0;
        plotting_options.plot_coil_current_weights = 0;
        plotting_options.plot_efields_coil = 0;
        plotting_options.plot_bfields_coil = 0;
        plotting_options.plot_coilSNR = 0;
        plotting_options.movie_coil_current = 0;
        plotting_options.plot_coil_number_flag = 0;
        plotting_options.plot_coil_setup_flag = 0;
        
        dgf_sphere_plot_snr(ult_output_values,path_options,plotting_options,snr_radius);
        
        plotting_options.plot_coilSNR = plot_coilSNR;
        plotting_options.plot_coil_current_patterns = plot_coil_current_patterns;
        plotting_options.plot_coil_current_weights = plot_coil_current_weights;
        plotting_options.plot_efields_coil = plot_efields_coil;
        plotting_options.plot_bfields_coil = plot_bfields_coil;
        plotting_options.plot_coilSNR = plot_coilSNR;
        plotting_options.movie_coil_current = movie_coil_current;
        plotting_options.plot_coil_number_flag = plot_coil_number_flag;
        plotting_options.plot_coil_setup_flag = plot_coil_setup_flag;
        
    end
    
    if plot_coilSNR && (~compute_coilSNR_flag)
        [coilSNRfile,coilSNRpath] = uigetfile([ circcoildir 'coils_snr_results/*.mat'],'Select the MAT-file with coilSNR results');
        load(fullfile(coilSNRpath,coilSNRfile));
        
        coil_output_values = struct(...
            'g_opts',g_opts,...
            's_opts',s_opts,...
            'nf',nf,...
            'np',np,...
            'snr_matrix_size',[2*tempindex_r 2*tempindex_c],...
            'whichcurrents',whichcurrents,...
            'fieldstrength',fieldstrength,...
            'acceleration',acceleration,...
            'snr_ult',[],...
            'g_ult',[],...
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
            'weights',[],...
            'weights_coil',weights_coil,...
            'currentpatternmatrixsize',currentpatternmatrixsize,...
            'currentpattern',[],...
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
            'efield_ult_set_fullfov',[],...
            'bfield_ult_set_fullfov',[],...
            'efield_coil_set_fullfov',efield_coil_set_fullfov,...
            'bfield_coil_set_fullfov',bfield_coil_set_fullfov,...
            'PRdivPL',NaN);
        
        if plot_regions_fov
            mu = 4*pi*1e-7;         % permeability of free space [Wb][A^-1][m^-1]
            c = 3e8;                % speed of light [m][s]
            epsilon_0 = 1/(mu*c^2); % permittivity [C][V^-1] = [s^4][A^2][m^-2][kg^-2]
            
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
            title(['Voxel Resolution = ' sprintf('%0.2f',1000*fovf/nf) ' [mm]'],'FontSize',16)
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
        
        plotting_options.plot_ultSNR = 0;
        plotting_options.plot_ult_current_patterns = 0;
        plotting_options.plot_ult_current_weights = 0;
        plotting_options.plot_efields_ult = 0;
        plotting_options.plot_bfields_ult = 0;
        plotting_options.plot_ultSNR = 0;
        plotting_options.movie_ult_current = 0;
        
        dgf_sphere_plot_snr(coil_output_values,path_options,plotting_options,snr_radius);
        
        plotting_options.plot_ultSNR = plot_ultSNR;
        plotting_options.plot_ult_current_patterns = plot_ult_current_patterns;
        plotting_options.plot_ult_current_weights = plot_ult_current_weights;
        plotting_options.plot_efields_ult = plot_efields_ult;
        plotting_options.plot_bfields_ult = plot_bfields_ult;
        plotting_options.plot_ultSNR = plot_ultSNR;
        plotting_options.movie_ult_current = movie_ult_current;
        
    end
    
    if save_log_file
        diary off
    end
end
% end

% END DGF_sphere_snr_batch.m
