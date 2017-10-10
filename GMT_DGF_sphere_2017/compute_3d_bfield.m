function [b1_basis,b1_coils] = compute_3d_bfield(x_fov,y_fov,z_fov,mask_reg1,mask_reg2,mask_reg3,mask_reg4,s_opts)

current_type = s_opts.current_type;
matrix_size = s_opts.matrix_size;
fieldstrength = s_opts.fieldstrength;
lmax = s_opts.lmax;
currentradius = s_opts.currentradius;
radius3 = s_opts.radius3;
radius2 = s_opts.radius2;
radius1 = s_opts.radius1;
ncoils = s_opts.ncoils;
coil_rotations = s_opts.coil_rotations;
coil_radii = s_opts.coil_radii;
coil_offsets = s_opts.coil_offsets;

epsilon_r4 = s_opts.epsilon_r4;
sigma_r4 = s_opts.sigma_r4;
epsilon_r3 = s_opts.epsilon_r3;
sigma_r3 = s_opts.sigma_r3;
epsilon_r2 = s_opts.epsilon_r2;
sigma_r2 = s_opts.sigma_r2;

mu = 4*pi*1e-7;         % permeability of free space [Wb][A^-1][m^-1]
c = 3e8;                % speed of light [m][s]
epsilon_0 = 1/(mu*c^2); % permittivity [C][V^-1] = [s^4][A^2][m^-2][kg^-2]

omega = 2*pi*42.576e6*fieldstrength; % Larmor frequency [Hz]

epsilon_r4 = epsilon_r4*epsilon_0;
epsilon_r3 = epsilon_r3*epsilon_0;
epsilon_r2 = epsilon_r2*epsilon_0;

k_0_squared = omega*omega*epsilon_0*mu;
k_0 = sqrt(k_0_squared);

k_4_squared = omega*mu*(omega*epsilon_r4+1i*sigma_r4);
k_4 = sqrt(k_4_squared);

k_3_squared = omega*mu*(omega*epsilon_r3+1i*sigma_r3);
k_3 = sqrt(k_3_squared);

k_2_squared = omega*mu*(omega*epsilon_r2+1i*sigma_r2);
k_2 = sqrt(k_2_squared);

if s_opts.compute_for_finite_arrays
    drfac = pi/180;
    rot_coil = coil_rotations*drfac;
end

% mask_reg1 = mask_reg1(:).';
mask_reg2 = mask_reg2(:).';
mask_reg3 = mask_reg3(:).';
mask_reg4 = mask_reg4(:).';

if length(current_type) > 1
    numbasis = 2*((lmax + 1)^2 - 1);
    counterincrement = 1;
else
    numbasis = (lmax + 1)^2 - 1;
    counterincrement = 0;
end
    
xset = x_fov(:).'; % vector of size [# of voxels, 1]
yset = y_fov(:).';
zset = z_fov(:).';
nvoxels = length(xset);

% -- constant factors depending on radius -- %
k_0_rad_b = k_0*currentradius;
besselnorm_rad_b = sqrt(pi/(2*k_0_rad_b));

% ele_scaling = omega*mu*k_0*currentradius*currentradius; % scaling factor for the Electric field (RL fixed on 6/11/2012)

mag_scaling_r4 = -1i*mu*k_4*k_0*currentradius*currentradius; % scaling factor for the Magnetic field (RL fixed on 6/11/2012)
mag_scaling_r3 = -1i*mu*k_3*k_0*currentradius*currentradius; 
mag_scaling_r2 = -1i*mu*k_2*k_0*currentradius*currentradius; 
mag_scaling_r1 = -1i*mu*k_0*k_0*currentradius*currentradius;

% Convert in spherical coordinates
rset = sqrt(xset.^2 + yset.^2 + zset.^2);   % rho
costhetaset = (zset./rset);                   % cos[theta]
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

krset = zeros(size(rset));
for irho = 1:length(rset)
    if rset(irho) < radius3
        krset(irho) = k_4*rset(irho);
    elseif (rset(irho) < radius2) && (rset(irho) >= radius3)
        krset(irho) = k_3*rset(irho);
    elseif (rset(irho) < radius1) && (rset(irho) >= radius2)
        krset(irho) = k_2*rset(irho);
    else
        krset(irho) = k_0*rset(irho);
    end
end
besselnorm = sqrt(pi./(2.*krset));

% initialize output variables
b1_allmodes = zeros(numbasis,nvoxels,3);

counter = 1;
for l = 1:lmax, % for l=0 the T matrix is empty and there would be no contributions
    %             disp(['      l = ' num2str(l)])
    lnorm = sqrt(l*(l+1));
    legendrenorm = sqrt((2*l + 1)/(4*pi));
    legendrenorm_minus1 = sqrt((2*l - 1)/(4*pi));
    legendrefunctions = legendre(l,costhetaset,'sch'); % # row is m = 0,...l ; # col is R
    if l>0,
        legendrefunctions_lminus1 = (legendre(l-1,costhetaset,'sch'));
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
        compute_reflection_and_transmission_coef(l,k_0,k_2,radius1);
    
    [R_H_P_2, R_H_F_2, R_V_P_2, R_V_F_2, T_H_P_2, T_H_F_2, T_V_P_2, T_V_F_2] = ...
        compute_reflection_and_transmission_coef(l,k_2,k_3,radius2);
    
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
        
        T = ((-1)^(1-m))*Ttemp; % RL Sept 25, 2017: I don't remember the reason for this
        
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
        
        if rset < eps % resolve singularity near the origin
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
        
        T_hat = T;
        if current_type == 1,
            T_hat = [1 0]*T;
        end
        if current_type == 2,
            T_hat = [0 1]*T;
        end
        
        % REGION 4
%         b1_allmodes(counter:counter+counterincrement,mask_reg4,1) = mag_scaling_r4*(T_hat*(repmat([D_m_41; D_n_41],1,nvoxels).*[N_x(mask_reg4); M_x(mask_reg4)]));
%         b1_allmodes(counter:counter+counterincrement,mask_reg4,2) = mag_scaling_r4*(T_hat*(repmat([D_m_41; D_n_41],1,nvoxels).*[N_y(mask_reg4); M_y(mask_reg4)]));
%         b1_allmodes(counter:counter+counterincrement,mask_reg4,3) = mag_scaling_r4*(T_hat*(repmat([D_m_41; D_n_41],1,nvoxels).*[N_z(mask_reg4); M_z(mask_reg4)]));
        b1_allmodes(counter:counter+counterincrement,mask_reg4,1) = mag_scaling_r4*(T_hat*[D_m_41*N_x(mask_reg4); D_n_41*M_x(mask_reg4)]);
        b1_allmodes(counter:counter+counterincrement,mask_reg4,2) = mag_scaling_r4*(T_hat*[D_m_41*N_y(mask_reg4); D_n_41*M_y(mask_reg4)]);
        b1_allmodes(counter:counter+counterincrement,mask_reg4,3) = mag_scaling_r4*(T_hat*[D_m_41*N_z(mask_reg4); D_n_41*M_z(mask_reg4)]);
        
        % REGION 3
        b1_allmodes(counter:counter+counterincrement,mask_reg3,1) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_x(mask_reg3) + D_m_31*N_x(mask_reg3)); (B_n_31*M3_x(mask_reg3) + D_n_31*M_x(mask_reg3)) ]);
        b1_allmodes(counter:counter+counterincrement,mask_reg3,2) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_y(mask_reg3) + D_m_31*N_y(mask_reg3)); (B_n_31*M3_y(mask_reg3) + D_n_31*M_y(mask_reg3)) ]);
        b1_allmodes(counter:counter+counterincrement,mask_reg3,3) = mag_scaling_r3*(T_hat*[ (B_m_31*N3_z(mask_reg3) + D_m_31*N_z(mask_reg3)); (B_n_31*M3_z(mask_reg3) + D_n_31*M_z(mask_reg3)) ]);
        
        % REGION 2
        b1_allmodes(counter:counter+counterincrement,mask_reg2,1) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_x(mask_reg2) + D_m_21*N_x(mask_reg2)); (B_n_21*M3_x(mask_reg2) + D_n_21*M_x(mask_reg2)) ]);
        b1_allmodes(counter:counter+counterincrement,mask_reg2,2) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_y(mask_reg2) + D_m_21*N_y(mask_reg2)); (B_n_21*M3_y(mask_reg2) + D_n_21*M_y(mask_reg2)) ]);
        b1_allmodes(counter:counter+counterincrement,mask_reg2,3) = mag_scaling_r2*(T_hat*[ (B_m_21*N3_z(mask_reg2) + D_m_21*N_z(mask_reg2)); (B_n_21*M3_z(mask_reg2) + D_n_21*M_z(mask_reg2)) ]);
        
        counter =  counter + 1 + counterincrement;
%         disp(['counter = ' num2str(counter)])
    end
end


% B1 for basis functions
if s_opts.compute_for_basis_functions
    if length(current_type) > 1
        b1_basis = zeros(2*s_opts.max_basis_functions,nvoxels,3);
        b1_basis = b1_allmodes(1:2*s_opts.max_basis_functions,:,:);
    else
        b1_basis = zeros(s_opts.max_basis_functions,nvoxels,3);
        b1_basis = b1_allmodes(1:s_opts.max_basis_functions,:,:);
    end
else
    b1_basis = [];
end

% B1 for finite arrays
if s_opts.compute_for_finite_arrays
    
    b1_coils = zeros(ncoils,nvoxels,3);
    if length(current_type) > 1
        b1_mag_modes = b1_allmodes(1:2:end,:,:);
        mag_modes = numbasis/2;
    else
        mag_modes = numbasis;
    end
    
    W_coil = zeros([((lmax + 1)^2 - 1) ncoils]);
    
    % define quantities constant for all coils
    costheta_coil_z = coil_offsets(1)/currentradius;
    theta_coil_z = acos(costheta_coil_z);
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
        
        if s_opts.reverse_m_flag
            W_coil_i = reverse_m_index(W_coil_i,lmax,current_type,1); % RL 25 June 2017: now it is a function
        end
        W_coil(:,icoil) = W_coil_i; % RL 16 Oct 2014 (THIS IS CORRECT, BUT CAUSES A PROBLEM WITH THE DISPLAY OF COIL CURRENTS)
        
        W_coil_i = W_coil_i(:);
                
        b1_coils(icoil,:,:) = reshape((W_coil_i.')*reshape(b1_mag_modes,[mag_modes nvoxels*3]),[nvoxels 3]);
    end
else
    b1_coils = [];
end

    
    
    


