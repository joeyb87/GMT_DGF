% -------------------------------------------------------------------------
%   DGF_sphere_gmt_batch.m
%
%   Global Maxwell Tomography (GMT) for a multi-layer spherical phantom
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

disp('DGF_sphere_gmt_batch');
disp(datetime)
disp('    --> Initialization...');
tic

% -------------------------------------------------------------------------
%  SET COMMONLY USED SIMULATION LOOP PARAMETERS
% -------------------------------------------------------------------------

%--- frequencies and fields ---%
% whichfrequencies = [1 32 64 96 128 192 256 300 350 400]*1e6;
whichfrequencies = [128]*1e6;
fieldstrength_set = whichfrequencies/42.576e6;
fieldstrength_set = [7]; %[1 3 5 7 9 11];
nfields = length(fieldstrength_set);

matrix_size = [32 32 32]; % [x y z]

%--- FOUR-LAYER GEOMETRY AND TISSUE PROPERTIES SETTINGS ---%

%__________________________________________________________%
% ** INNER MOST LAYER - REGION 4 (0 < radius < radius3) ** %
%----------------------------------------------------------%

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

%__________________________________________________________%
% **  REGION 3 (radius3 < radius < radius2) ** %
%----------------------------------------------------------%

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

%__________________________________________________________%
% **  REGION 2 (radius2 < radius < radius1) ** %
%----------------------------------------------------------%

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

%--- allowed current types ---% (1 = magnetic-dipole, 2 = electric-dipole)
% current_type =      {1   1   [1 2]};
current_type =      [1 2];
% current_type =      {2};

%--- order of the expansion ---%
lmax = 10; % # of modes = 2(lmax + 1)^2
% NOTE: lmax must be sufficiently large to define coils, so use "max_basis_functions" to
%       restrict the number of modes used and avoid saving very large files

% -------------------------------------------------------------------------
%  SET SIMULATION FLAGS AND OPTIONS
% -------------------------------------------------------------------------

measurement_type = 'b1plus';
% 'b1plus' --> use only B1+ (with relative phases) in the cost function
% 'mrsignal' --> use only the MR signal in the cost function
% 'fullmeasurement' --> use both B1+ (with relative phases) and the MR signal in the cost function

shimming_type = 'none';
% 'none' --> transmit elements are used individually and driven with unit weights
% 'all_but_one' --> virtual excitations are created in an all-but-one fashion enforcing constructive interference at the center [DEFAULT]

initial_guess = 'homogeneous';
% 'homogeneous' --> average EP value at every voxel

%--- general ---%
warningsoffflag = 1;                    % 1 --> set off some MATLAB warnings
reverse_m_flag = 1;                     % 1 --> reverse mode order with respect to m

%--- simulation options ---%
compute_for_basis_functions = 1;        % 1 --> perform GMT using current modes as hypothetical elements of a transmit/receive array
max_basis_functions = 16;               % defines how many basis functions are used
compute_for_finite_arrays = 1;          % 1 --> perform GMT with a finite transmit/receive array
compute_current_patterns = 0;           % 1 --> compute current patterns with same combination weights used for the transmit fields
add_noise = 1;                          % 1 --> add white gaussian noise to the synthetic measurements
peak_SNR = 150;                         % establish the peak SNR, which defines how much noise will be added to the synthetic measurements
preset_coil_geometry_flag = 0;          % 1 --> coil geometry is loaded from file
preset_coil_weights_flag = 0;           % 1 --> coil weights are loaded from a stored MAT file (will bypass shimming_type)
transmit_subset = {'all'};             % defines which coils are going to be used for transmission
receive_subset =  {'all'};              % defines which coils are going to be used for reception

%--- saved variables ---%
save_currents_flag = 1;                 % 1 --> save current patterns
save_basis_field_patterns_flag = 0;     % 1 --> save modes' raw B field patterns
save_coil_field_patterns_flag = 0;      % 1 --> save coils' raw B field patterns
save_individual_basis_b1plus = 0;       % 1 --> save shimmed B1+ for each basis function used (these are synthetic measurements used in GMT)
save_individual_basis_signal = 0;       % 1 --> save MR signal for each basis function used (these are synthetic measurements used in GMT)
save_individual_coils_b1plus = 0;       % 1 --> save shimmed B1+ for each coil in the array (these are synthetic measurements used in GMT)
save_individual_coils_signal = 0;       % 1 --> save MR signal for each coil in the array (these are synthetic measurements used in GMT)
save_basis_weights = 0;                 % 1 --> save the combination weights used for the modes' transmit fields
save_coil_weights = 0;                  % 1 --> save the combination weights used for the coils' transmit fields

%--- plotting ---%
plot_stored_results = 0;                % 1 --> plot results loaded from MAT file (using plot options specified here and not in the stored result file; WARNING: no computation performed if this flag = 1)
plot_gmt_all = 1;                       % 1 --> plot all GMT results
plot_ep_maps = 0;                       % 1 --> plot reconstructed maps of electrical properties
plot_error = 1;                         % 1 --> 2x2 plot with initial guess, true properties, reconstructed properties and signed error
plot_currents = 0;                      % 1 --> plot current patterns computed by combining the modes with the same combination weights used for the transmit fields
plot_currents_3D = 1;                   % 1 --> display currents pattern on tridimensional surfaces (only if current patterns are plotted)
plot_net_field_patterns = 0;            % 1 --> plot net B1+ field patterns (abs and phase of net unshimmed field) for both basis and array
plot_individual_coils_fields = 0;       % 1 --> plot B1+ field patterns for each coil in the array (abs and phase of unshimmed coil field)
plot_excitation_patterns = 0;           % 1 --> plot B1+ field patterns (both basis and array) for each excitation used in GMT (e.g., all but one combinations)
plot_object_setup_flag = 0;             % 1 --> figure with layers and current surface
plot_coil_setup = 0;                    % 1 --> figure with sphere and array
plot_coil_numbers = 0;                  % 1 --> add coil number next to array elements

%--- general plotting options ---%
user_label = 'RL';                      % user signature that will be added on some of the plots
plot_font_size = 18;                    % FontSize used for some of the plots
allow_windowing_flag = 0;               % 1 --> shows sliding bars to adjust windowing on figure windows (i.e. call to imgui, which can be undone for the active figure by entering imgui in the command window)
save_plots = 0;                         % ## TO BE DONE ## 1 --> save all figures generated in separate files using export_fig
save_plots_format = 'png';              % available formats: ai, bmp, emf, eps, fig, jpg, pdf, png, tiff
save_plots_resolution = 150;            % specify the resolution of the exported figure

%--- current patterns plotting options ---%
real_part_flag = 1;                     % 1 --> plot the real part of the complex-valued current patterns (0 --> plot imaginary part)
show_axes = 1;                          % 1 --> show Cartesian axes in the 3D plot
axes_color = [0 0 0];                   % 1 --> color of Cartesian axes in the 3D plot
plot_bkgcolor = [0.5 0.5 0.5];          % background color for 3D plots [1 1 1] = white, [0.5 0.5 0.5] = gray, [0 0 0] = black
currentpatternmatrixsize = [52 52];     % number of points at which the current vectors are calculated
imagebackgroundflag = 1;                % ## TO BE DONE ## 1 --> plot current patterns with a background image showing current magnitude
patternphase = 90;                      % ## TO BE DONE ## phase of current or field pattern plotted (degrees)
currentarrowswidth = 1.1;               % "LineWidth" for quiver plot of current patterns
currentarrowscolor = 'k';               % "Color" for quiver plot of current patterns
spherefacecolor = [0.89 0.89 0.89];     % "FaceColor" color of the sphere in 3D current plots. [0.89 0.89 0.89] = Black&White figures, [1 0.46 0] = Ultimate basis set, [0.9 1 1] = Light blue sphere
sphereedgecolor = 'none';               % "EdgeColor" color of the patches' borders in 3D current plots
spherefacealpha = 1;                    % "FaceAlpha" level of transparency of the patches in 3D current plots
sphereview = [-37.5-180 30];            % azimuth and elevation of the viewpoint in 3D current plots

%--- EM fields plotting options ---%
scale_factor = 2;                       % new matrix size will be scale_factor*matrix_size
interp_method = 'bilinear';             % interpolation method for scaling the original image (nearest, bilinear, bicubic)
display_range_ult = [0 1.4e-5];         % display range of imshow (for coil fields), specified as a two-element vector [LOW HIGH]; [] --> default
display_range_coil = [0 1.4e-5];        % display range of imshow (for coil fields), specified as a two-element vector [LOW HIGH]; [] --> default
plot_profiles = 1;                      % plot vertical and horizontal profiles
logplot = 1;                            % plot absolute value of fields in log10 scale

% -------------------------------------------------------------------------
% DEFINE DIRECTORIES PATHS
% -------------------------------------------------------------------------

home_dir     = './';
commondir   = [home_dir 'common_routines/'];
plotdir     = [home_dir 'plot_routines/'];
resultsdir  = [home_dir 'results/'];
geometrydir = [home_dir 'coil_packing/'];
tissuedir   = [home_dir 'tissue_properties/'];

addpath(plotdir);
addpath(commondir);

path_options = struct(...
    'home_dir',home_dir,...
    'commondir',commondir,...
    'plotdir',plotdir,...
    'resultsdir',resultsdir,...
    'geometrydir',geometrydir,...
    'tissuedir',tissuedir);

if plot_gmt_all
    plot_ep_maps = 1;
    plot_error = 1;
    plot_currents = 1;
    plot_currents_3D = 1;
    plot_ult_field_patterns = 1;
    plot_array_field_patterns = 1;
    plot_individual_coils_fields = 1;
    plot_object_setup_flag = 1;
    plot_coil_setup = 1;
end

plotting_options = struct(...
    'plot_ep_maps',plot_ep_maps,...
    'plot_error',plot_error,...
    'plot_currents',plot_currents,...
    'plot_currents_3D',plot_currents_3D,...
    'plot_net_field_patterns',plot_net_field_patterns,...
    'plot_individual_coils_fields',plot_individual_coils_fields,...
    'plot_excitation_patterns',plot_excitation_patterns,...
    'plot_object_setup_flag',plot_object_setup_flag,...
    'plot_coil_setup',plot_coil_setup,...
    'plot_coil_numbers',plot_coil_numbers,...
    'user_label',user_label,...
    'plot_font_size',plot_font_size,...
    'allow_windowing_flag',allow_windowing_flag,...
    'save_plots',save_plots,...
    'save_plots_format',save_plots_format,...
    'save_plots_resolution',save_plots_resolution',...
    'real_part_flag',real_part_flag,...
    'show_axes',show_axes,...
    'axes_color',axes_color,...
    'plot_bkgcolor',plot_bkgcolor,...
    'currentpatternmatrixsize',currentpatternmatrixsize,...
    'imagebackgroundflag',imagebackgroundflag,...
    'patternphase',patternphase,...
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
    'plot_profiles',plot_profiles,...
    'logplot',logplot);

results_file = [];
if plot_stored_results % Plot stored results, no computation
    
    disp('** WARNING: plotting stored results, no computation performed **')
    compute_for_basis_functions = 0;
    compute_for_finite_arrays = 0;
    compute_current_patterns = 0;
    
    if results_file
        gmt_results = load(results_file,gmt_results)
    else
        [results_file,results_path] = uigetfile([ basissetdir 'ultimate_snr_results/*.mat'],'Select the MAT-file with ultSNR results');
        gmt_results = load(fullfile(results_path,results_file),gmt_results);
    end
    dgf_sphere_gmt_plot(gmt_results,plotting_options,path_options)
    
else
    
    if (compute_for_basis_functions == 0) && (compute_for_finite_arrays == 0)
        disp('-------------------------------------------------------');
        disp('## ERROR: set at least one computation flag   ');
        disp('-------------------------------------------------------');
        keyboard
    end

    % -------------------------------------------------------------------------
    %  DEFINE ARRAY GEOMETRY AND LOAD SHIMMING COEFFICIENTS
    % -------------------------------------------------------------------------
    
    % coil_file = [];
    coil_file = [geometrydir '8_coils_mv.mat'];
    coil_weights_file = [];  % will be used if preset_coil_weights_flag = 1
    
    if compute_for_finite_arrays
        if preset_coil_geometry_flag
            if coil_file
                load(coil_file)
            else
                [coil_file,coil_path] = uigetfile('*.mat','Select the MAT-file with coil array geometry',geometrydir);
                load(fullfile(coil_path,coil_file));
            end
            ncoils = size(coil_rotations,1);
            coil_rotations = 180*coil_rotations/pi; % conversion from rad to deg
        else
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
            coil_radius_mul = 0.36;
            coil_surface_radii = ones(ncoils,1)*currentradius;
            coil_radii = coil_surface_radii*coil_radius_mul;
            coil_offsets = sqrt(coil_surface_radii.^2 - coil_radii.^2);
        end
        %--- check if all varaiables are there ---%
        if any(isnan(ncoils)) || isempty(coil_rotations) || any(isnan(coil_radius_mul))
            disp('-------------------------------------------------------');
            disp('**ERROR** unspecified circular surface coils variable');
            disp('-------------------------------------------------------');
            pause
        end
        if plot_coil_setup
            figure;
            plot_geometry_sphere(radius1,coil_radii,coil_offsets,coil_rotations,[0.1],[0],[0.1],3,'-',[1 0.46 0],0,plot_coil_numbers)
            plotting_options.plot_coil_setup = 0; % to avoid plotting it again
        end
        if preset_coil_weights_flag
            % the following will load the variable "coil_weights" with the shimming weights for the array elements.
            % if ~preset_coil_weights_flag then "b1_shim_coils" will be created after calculating the fields
            if coil_weights_file
                b1_shim_coils = load(coil_weights_file,coil_weights)
            else
                [coil_weights_file,weights_path] = uigetfile('*.mat','Select the MAT-file with coil combination weights',geometrydir);
                b1_shim_coils = load(fullfile(weights_path,coil_weights_file),coil_weights);
            end
            shimming_type = 'all_but_one'; % Setting it to the default value, to be used eventually only for basis functions
        else
            b1_shim_coils = [];
        end
    else
        ncoils = NaN;
        coil_rotations = [];
        coil_radii = [];
        coil_offsets = [];
        b1_shim_coils = [];
    end
    b1_shim_basis = []; % It will be updated after calculating the fields
    
    % -------------------------------------------------------------------------
    %  GROUPING SIMULATION SETTINGS AND FLAGS
    % -------------------------------------------------------------------------
    
    if save_basis_field_patterns_flag && ~(compute_for_basis_functions)
        save_basis_field_patterns_flag = 0;
        disp('** WARNING: cannot save modes'' net E and B fields because they are not computed **')
    end
    if save_coil_field_patterns_flag && ~(compute_for_finite_arrays)
        save_coil_field_patterns_flag = 0;
        disp('** WARNING: cannot save array''s net E and B fields because they are not computed **')
    end
    if save_individual_coils_b1plus && ~(compute_for_finite_arrays)
        save_individual_coils_b1plus = 0;
        disp('** WARNING: cannot save coils'' individual B1+ fields because they are not computed **')
    end
    if save_individual_coils_signal && ~(compute_for_finite_arrays)
        save_individual_coils_signal = 0;
        disp('** WARNING: cannot save coils'' individual signals because they are not computed **')
    end
    if save_basis_weights && ~(compute_for_basis_functions)
        save_basis_weights = 0;
        disp('** WARNING: cannot save modes'' weights because they are not computed **')
    end
    if save_coil_weights && ~(compute_for_finite_arrays)
        save_coil_weights = 0;
        disp('** WARNING: cannot save arrays'' weights because they are not computed **')
    end
    
    if save_currents_flag && ~(compute_current_patterns)
        save_currents_flag = 0;
        disp('** WARNING: cannot save current patterns because they are not computed **')
    end
    
    simulation_options = struct(...
        'reverse_m_flag',reverse_m_flag,...
        'matrix_size',matrix_size,...
        'radius3',radius3,...
        'tissue_r4',tissue_r4,...
        'epsilon_r4',epsilon_r4,...
        'sigma_r4',sigma_r4,...
        'radius2',radius2,...
        'tissue_r3',tissue_r3,...
        'epsilon_r3',epsilon_r3,...
        'sigma_r3',sigma_r3,...
        'radius1',radius1,...
        'tissue_r2',tissue_r2,...
        'epsilon_r2',epsilon_r2,...
        'sigma_r2',sigma_r2,...
        'currentradius',currentradius,...
        'current_type',current_type,...
        'lmax',lmax,...
        'measurement_type',measurement_type,...
        'shimming_type',shimming_type,...
        'initial_guess',initial_guess,...
        'compute_for_basis_functions',compute_for_basis_functions,...
        'max_basis_functions',max_basis_functions,...
        'compute_for_finite_arrays',compute_for_finite_arrays,...
        'compute_current_patterns',compute_current_patterns,...
        'add_noise',add_noise,...
        'peak_SNR',peak_SNR,...
        'save_currents_flag',save_currents_flag,...
        'save_basis_field_patterns_flag',save_basis_field_patterns_flag,...
        'save_coil_field_patterns_flag',save_coil_field_patterns_flag,...
        'save_individual_basis_b1plus',save_individual_coils_b1plus,...
        'save_individual_basis_signal',save_individual_coils_signal,...
        'save_individual_coils_b1plus',save_individual_coils_b1plus,...
        'save_individual_coils_signal',save_individual_coils_signal,...
        'save_basis_weights',save_basis_weights,...
        'save_coil_weights',save_coil_weights,...
        'transmit_subset',transmit_subset,...
        'receive_subset',receive_subset,...
        'ncoils',ncoils,...
        'coil_rotations',coil_rotations,...
        'coil_radii',coil_radii,...
        'coil_offsets',coil_offsets,...
        'b1_shim_coils',b1_shim_coils,...
        'b1_shim_basis',b1_shim_basis);
    
    % -------------------------------------------------------------------------
    %  BEGIN SIMULATION
    % -------------------------------------------------------------------------
    
    if warningsoffflag,
        warning('off','MATLAB:nearlySingularMatrix')
        warning('off','MATLAB:divideByZero')
    end
    
    if max_basis_functions > 128
        disp('** WARNING: reduce the number of basis functions to avoid saving large files **')
    end
    
    % Error handling
    dbstop if error
    
    fieldcount = 1;
    timecheck = toc;
    disp(['Initialization completed (', num2str(floor(timecheck/3600)), ' h ' num2str(floor(mod(timecheck,3600)/60)) ' m ' num2str(round(mod(mod(timecheck,3600),60)),1) ' s)']);
    disp('    --> Starting GMT simulation for...');
    for ifield = 1:nfields
        fieldstrength = fieldstrength_set(ifield);
        simulation_options.fieldstrength = fieldstrength; % added to the structure
        disp(['        B0 = ',num2str(fieldstrength),' Tesla']);
        
        if isnan(sigma_r4) || isnan(epsilon_r4),
            tissuefile_r4 = [tissue_r4 '.mat'];
            load(tissuefile_r4);
            fieldset = tissueproperties(1,:)./42.576e6;
            sigma_set = tissueproperties(2,:);
            epsilon_rel_set = tissueproperties(3,:);
            epsilon_r4 = spline(fieldset,epsilon_rel_set,fieldstrength);
            sigma_r4 = spline(fieldset,sigma_set,fieldstrength);
            clear tissueproperties
            simulation_options.epsilon_r4 = epsilon_r4;
            simulation_options.sigma_r4 = sigma_r4;
        end
        if isnan(sigma_r3) || isnan(epsilon_r3),
            tissuefile_r3 = [tissue_r3 '.mat'];
            load(tissuefile_r3);
            fieldset = tissueproperties(1,:)./42.576e6;
            sigma_set = tissueproperties(2,:);
            epsilon_rel_set = tissueproperties(3,:);
            epsilon_r3 = spline(fieldset,epsilon_rel_set,fieldstrength);
            sigma_r3 = spline(fieldset,sigma_set,fieldstrength);
            clear tissueproperties
            simulation_options.epsilon_r3 = epsilon_r3;
            simulation_options.sigma_r3 = sigma_r3;
        end
        if isnan(sigma_r2) || isnan(epsilon_r2),
            tissuefile_r2 = [tissue_r2 '.mat'];
            load(tissuefile_r2);
            fieldset = tissueproperties(1,:)./42.576e6;
            sigma_set = tissueproperties(2,:);
            epsilon_rel_set = tissueproperties(3,:);
            epsilon_r2 = spline(fieldset,epsilon_rel_set,fieldstrength);
            sigma_r2 = spline(fieldset,sigma_set,fieldstrength);
            clear tissueproperties
            simulation_options.epsilon_r2 = epsilon_r2;
            simulation_options.sigma_r2 = sigma_r2;
        end
        
        %% create object masks and voxel coordinates
        if fieldcount == 1,
            [x_fov,y_fov,z_fov,mask_reg1,mask_reg2,mask_reg3,mask_reg4] = create_3d_fov_sphere(simulation_options);
            if plot_object_setup_flag
                plot_object_setup(simulation_options)
            end
        end
        
        %% compute B fields for modes/coils
        [bfield_basis,bfield_coils] = compute_3d_bfield(x_fov,y_fov,z_fov,mask_reg1,mask_reg2,mask_reg3,mask_reg4,simulation_options);
        
        % compute B1+ fields for modes/coils
        if ~isempty(bfield_basis)
            b1plus_basis = bfield_basis(:,:,1) + 1i*bfield_basis(:,:,2);
        else
            b1plus_basis = [];
        end
        if ~isempty(bfield_coils)
            b1plus_coils = bfield_coils(:,:,1) + 1i*bfield_coils(:,:,2);
        else
            b1plus_coils = [];
        end
        
        %% --- TEMPORARY (start) --- %
        % CHECK B1+ FIELDS
        check_fields_flag = 1;
        if check_fields_flag
            if ~isempty(b1plus_basis)
                netb1p_basis = squeeze(sum(b1plus_basis,1));
                netb1p_basis = reshape(netb1p_basis,matrix_size);
                view_3d_field(abs(netb1p_basis));
            end
            if ~isempty(b1plus_coils)
                netb1p_coils = squeeze(sum(b1plus_coils,1));
                netb1p_coils = reshape(netb1p_coils,matrix_size);
                view_3d_field(abs(netb1p_coils));
            end
        end
        % --- TEMPORARY (end) --- %       
        
        %% compute b1_shim
        [b1shim_basis,b1shim_coils] = compute_b1shims(b1plus_basis,b1plus_coils,simulation_options);
        
        %% generate synthetic measurements
        [shimmed_b1plus_basis,shimmed_b1plus_coils,signal_basis,signal_coils] = generate_measurements(b1plus_basis,b1plus_coils,b1shim_basis,b1shim_coils,simulation_options);
        
        %% generate initial guess
        [eps_r_initial,sigma_initial] = generate_initial_guess(simulation_options);
        
        
        % ## DONE ### 1) create object masks and voxel coordinates
        % ## DONE ### 2) compute B fields for modes/coils
        % 3) compute b1_shim
        % 4) generate measurements
        % 5) generate initial guess
        % 6) run GMT reconstruction iteratively
        
        
        timecheck = toc;
        disp(['-----------------------------'])
        disp(['GMT simulation completed for B0 = ',num2str(fieldstrength),' Tesla (', num2str(floor(timecheck/3600)), ' h ' num2str(floor(mod(timecheck,3600)/60)) ' m ' num2str(round(mod(mod(timecheck,3600),60)),1) ' s)']);
        disp(['    # of Iterations: ',num2str(niterations)])
        disp(['-----------------------------'])
        fieldcount = fieldcount + 1;
    end
end

% END - DGF_sphere_gmt_batch.m
