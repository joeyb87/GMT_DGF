function   [Efieldpattern,Bfieldpattern,Efieldpattern_coil,Bfieldpattern_coil] = ...
    calculate_field_patterns_function(...
    whichcurrents,s_opts,T_matrix_set,l,m,currentweights,W_coil,currentweights_coil,test_1,test_2,counter_fields,...
    radius3,radius2,radius1,k_0_rad_b,mask_r4,mask_r3,mask_r2,mask_r1,...
    legendrefunctions,legendrefunctions_lminus1,...
    mag_scaling_r1,mag_scaling_r2,mag_scaling_r3,mag_scaling_r4,ele_scaling,...
    k_0_squared,k_0,k_2_squared,k_2,k_3_squared,k_3,k_4_squared,k_4,...
    fieldpatternmatrixsize,Efieldpattern,Bfieldpattern,Efieldpattern_coil,Bfieldpattern_coil,...
    fieldrho,cosfieldtheta,fieldtheta,fieldphi,sinfieldtheta,sinfieldphi,...
    cosfieldphi,cotfieldtheta,cscfieldtheta,cos2fieldtheta)

% allocate arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Efpt = zeros([fieldpatternmatrixsize 3]); % last dimension is x,y,z components
Bfpt = Efpt;
if s_opts.compute_coil_field_patterns
    Efpt_coil = Efpt;
    Bfpt_coil = Bfpt;
end
E_x = zeros([length(whichcurrents) fieldpatternmatrixsize]);
E_y = E_x;
E_z = E_x;
B_x = E_x;
B_y = E_x;
B_z = E_x;

rhat_x = sinfieldtheta.*cosfieldphi;
rhat_y = sinfieldtheta.*sinfieldphi;
rhat_z = cosfieldtheta;

krset = zeros(size(fieldrho));
krset(mask_r4) = k_4*fieldrho(mask_r4);
krset(mask_r3) = k_3*fieldrho(mask_r3);
krset(mask_r2) = k_2*fieldrho(mask_r2);
krset(mask_r1) = k_0*fieldrho(mask_r1);

besselnorm = sqrt(pi./(2.*krset));
besselnorm_rad_b = sqrt(pi/(2*k_0_rad_b));

lnorm = sqrt(l*(l+1));
legendrenorm = sqrt((2*l + 1)/(4*pi));
legendrenorm_minus1 = sqrt((2*l - 1)/(4*pi));
%     legendrefunctions = legendre(l,cosfieldtheta,'sch'); % # row is m = 0,...l ; # col is R
%     if l>0,
%         legendrefunctions_lminus1 = legendre(l-1,cosfieldtheta,'sch');
%         legendrefunctions_lminus1 = [legendrefunctions_lminus1; zeros(size(cosfieldtheta))];
%     end

% -- these two are for rho = sphere radius --
%             H_l_1_k_0_rad_a = besselnorm_rad_a*besselh(l+0.5,1,k_0_rad_a);
%             H_l_1_k_0_rad_a_prime = -k_0_rad_a*besselnorm_rad_a*besselh(l+1+0.5,1,k_0_rad_a) + (l + 1)*H_l_1_k_0_rad_a;
% -----------------------------

H_l_1_k_0_rad_b = besselnorm_rad_b*besselh(l+0.5,1,k_0_rad_b);
H_l_1_k_0_rad_b_prime = -k_0_rad_b*besselnorm_rad_b*besselh(l+1+0.5,1,k_0_rad_b) + (l + 1)*H_l_1_k_0_rad_b;

J_l_kr = besselnorm.*besselj(l+0.5,krset);
J_l_kr(fieldrho<eps) = (l==0);                            %** only j_0 survives at r=0
J_lplus1_kr = besselnorm.*besselj(l+1+0.5,krset);
J_lplus1_kr(fieldrho<eps) = 0;                           %** j_l+1 vanishes at r=0 for l>=0
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

D_m_41 = repmat(D_m_41,[fieldpatternmatrixsize 1]);
D_n_41 = repmat(D_n_41,[fieldpatternmatrixsize 1]);

B_m_31 = repmat(B_m_31,[fieldpatternmatrixsize 1]);
B_n_31 = repmat(B_n_31,[fieldpatternmatrixsize 1]);
D_m_31 = repmat(D_m_31,[fieldpatternmatrixsize 1]);
D_n_31 = repmat(D_n_31,[fieldpatternmatrixsize 1]);

B_m_21 = repmat(B_m_21,[fieldpatternmatrixsize 1]);
B_n_21 = repmat(B_n_21,[fieldpatternmatrixsize 1]);
D_m_21 = repmat(D_m_21,[fieldpatternmatrixsize 1]);
D_n_21 = repmat(D_n_21,[fieldpatternmatrixsize 1]);


lmul = sqrt((2*l + 1)*(l^2 - m^2)/(2*l - 1));

legdim = size(legendrefunctions);
if m > 0
    Y_l_m = ((-1)^m)*legendrenorm*(1/sqrt(2))*reshape(legendrefunctions((m+1),:,:,:),[legdim(2:end) 1]).*exp(1i*m*fieldphi);
elseif m == 0
    Y_l_m = legendrenorm*reshape(legendrefunctions(1,:,:,:),[legdim(2:end) 1]);
else % using Y_(l,-m) = (-1)^m*conj[Y_(l,m)]
    Y_l_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm*(1/sqrt(2))*reshape(legendrefunctions((abs(m)+1),:,:,:),[legdim(2:end) 1]).*exp(1i*abs(m)*fieldphi));
end
if m > 0
    Y_lminus1_m = ((-1)^m)*legendrenorm_minus1*(1/sqrt(2))*reshape(legendrefunctions_lminus1(m+1,:,:,:),[legdim(2:end) 1]).*exp(1i*m*fieldphi);
elseif m == 0
    Y_lminus1_m = legendrenorm_minus1*reshape(legendrefunctions_lminus1(1,:,:),[legdim(2:end) 1]);
else
    Y_lminus1_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm_minus1*(1/sqrt(2))*reshape(legendrefunctions_lminus1((abs(m)+1),:,:),[legdim(2:end) 1]).*exp(1i*abs(m)*fieldphi));
end

X_x = (1/lnorm).*((-m*cosfieldphi+1i*l*sinfieldphi).*cotfieldtheta.*Y_l_m ...
    - 1i*lmul*cscfieldtheta.*sinfieldphi.*Y_lminus1_m);
X_y = (1/lnorm).*((-m*sinfieldphi-1i*l*cosfieldphi).*cotfieldtheta.*Y_l_m ...
    + 1i*lmul*cscfieldtheta.*cosfieldphi.*Y_lminus1_m);
X_z = (1/lnorm).*m.*Y_l_m;

r_cross_X_x = (1/lnorm).*((m*sinfieldphi+1i*l*cos2fieldtheta.*cosfieldphi).*cscfieldtheta.*Y_l_m ...
    - 1i*lmul*cotfieldtheta.*cosfieldphi.*Y_lminus1_m);
r_cross_X_y = (1/lnorm).*((-m*cosfieldphi+1i*l*cos2fieldtheta.*sinfieldphi).*cscfieldtheta.*Y_l_m ...
    - 1i*lmul*cotfieldtheta.*sinfieldphi.*Y_lminus1_m);
r_cross_X_z = (-1i/lnorm).*(l*cotfieldtheta.*Y_l_m ...
    - lmul.*Y_lminus1_m);

if ndims(fieldtheta)>2 % therefore is not the transverse plane. RL: need to generalize this
    %     [smtheta_1, smtheta_2, smtheta_3] = find( (abs(fieldtheta) < 0.01) | ((abs(fieldtheta)-pi) < 0.01) );
    %     if m == 1
    %         X_x(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(cosfieldtheta(smtheta_1,smtheta_2,smtheta_3).*cosfieldphi(smtheta_1,smtheta_2,smtheta_3) - 1i*sinfieldphi(smtheta_1,smtheta_2,smtheta_3));
    %         X_y(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(cosfieldtheta(smtheta_1,smtheta_2,smtheta_3).*sinfieldphi(smtheta_1,smtheta_2,smtheta_3) + 1i*cosfieldphi(smtheta_1,smtheta_2,smtheta_3));
    %         X_z(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(-sinfieldtheta(smtheta_1,smtheta_2,smtheta_3));
    %         r_cross_X_x(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(-1i*cosfieldtheta(smtheta_1,smtheta_2,smtheta_3).*cosfieldphi(smtheta_1,smtheta_2,smtheta_3) - sinfieldphi(smtheta_1,smtheta_2,smtheta_3));
    %         r_cross_X_y(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(-1i*cosfieldtheta(smtheta_1,smtheta_2,smtheta_3).*sinfieldphi(smtheta_1,smtheta_2,smtheta_3) + cosfieldphi(smtheta_1,smtheta_2,smtheta_3));
    %         r_cross_X_z(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(+1i*sinfieldtheta(smtheta_1,smtheta_2,smtheta_3));
    %     elseif m == -1
    %         X_x(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(cosfieldtheta(smtheta_1,smtheta_2,smtheta_3).*cosfieldphi(smtheta_1,smtheta_2,smtheta_3) + 1i*sinfieldphi(smtheta_1,smtheta_2,smtheta_3));
    %         X_y(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(cosfieldtheta(smtheta_1,smtheta_2,smtheta_3).*sinfieldphi(smtheta_1,smtheta_2,smtheta_3) - 1i*cosfieldphi(smtheta_1,smtheta_2,smtheta_3));
    %         X_z(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(-sinfieldtheta(smtheta_1,smtheta_2,smtheta_3));
    %         r_cross_X_x(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(1i*cosfieldtheta(smtheta_1,smtheta_2,smtheta_3).*cosfieldphi(smtheta_1,smtheta_2,smtheta_3) - sinfieldphi(smtheta_1,smtheta_2,smtheta_3));
    %         r_cross_X_y(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(1i*cosfieldtheta(smtheta_1,smtheta_2,smtheta_3).*sinfieldphi(smtheta_1,smtheta_2,smtheta_3) + cosfieldphi(smtheta_1,smtheta_2,smtheta_3));
    %         r_cross_X_z(smtheta_1,smtheta_2,smtheta_3) = (1/2)*sqrt((2*l + 1)/(4*pi))*(-1i*sinfieldtheta(smtheta_1,smtheta_2,smtheta_3));
    %     else
    %         X_x(smtheta_1,smtheta_2,smtheta_3) = 0;
    %         X_y(smtheta_1,smtheta_2,smtheta_3) = 0;
    %         X_z(smtheta_1,smtheta_2,smtheta_3) = 0;
    %         r_cross_X_x(smtheta_1,smtheta_2,smtheta_3) = 0;
    %         r_cross_X_y(smtheta_1,smtheta_2,smtheta_3) = 0;
    %         r_cross_X_z(smtheta_1,smtheta_2,smtheta_3) = 0;
    %     end
    [smtheta] = find(abs(fieldtheta) < 0.01);
    [bgtheta] = find(squeeze(abs(abs(fieldtheta)-pi)) < 0.01);
    use_Hong_Hsi_derivation = 0; % 1 --> analytic derivation; 0 --> series expansion
    if use_Hong_Hsi_derivation
        if m == 1
            X_x(smtheta) = (1/2)*sqrt((2*l + 1)/(4*pi));
            X_y(smtheta) = 1i*(1/2)*sqrt((2*l + 1)/(4*pi));
            X_z(smtheta) = 0;
            r_cross_X_x(smtheta) = -1i*(1/2)*sqrt((2*l + 1)/(4*pi));
            r_cross_X_y(smtheta) = (1/2)*sqrt((2*l + 1)/(4*pi));
            r_cross_X_z(smtheta) = 0;
            
            X_x(bgtheta) = -(1/2)*sqrt((2*l + 1)/(4*pi))*(-1^(l+1));
            X_y(bgtheta) = -1i*(1/2)*sqrt((2*l + 1)/(4*pi))*(-1^(l+1));
            X_z(bgtheta) = 0;
            r_cross_X_x(bgtheta) = -1i*(1/2)*sqrt((2*l + 1)/(4*pi))*(-1^(l+1));
            r_cross_X_y(bgtheta) = (1/2)*sqrt((2*l + 1)/(4*pi))*(-1^(l+1));
            r_cross_X_z(bgtheta) = 0;
        elseif m == -1
            X_x(smtheta) = (1/2)*sqrt((2*l + 1)/(4*pi));
            X_y(smtheta) = -1i*(1/2)*sqrt((2*l + 1)/(4*pi));
            X_z(smtheta) = 0;
            r_cross_X_x(smtheta) = 1i*(1/2)*sqrt((2*l + 1)/(4*pi));
            r_cross_X_y(smtheta) = (1/2)*sqrt((2*l + 1)/(4*pi));
            r_cross_X_z(smtheta) = 0;
            
            X_x(bgtheta) = -(1/2)*sqrt((2*l + 1)/(4*pi))*(-1^(l+1));
            X_y(bgtheta) = 1i*(1/2)*sqrt((2*l + 1)/(4*pi))*(-1^(l+1));
            X_z(bgtheta) = 0;
            r_cross_X_x(bgtheta) = 1i*(1/2)*sqrt((2*l + 1)/(4*pi))*(-1^(l+1));
            r_cross_X_y(bgtheta) = (1/2)*sqrt((2*l + 1)/(4*pi))*(-1^(l+1));
            r_cross_X_z(bgtheta) = 0;
        else
            X_x(smtheta) = 0;
            X_y(smtheta) = 0;
            X_z(smtheta) = 0;
            r_cross_X_x(smtheta) = 0;
            r_cross_X_y(smtheta) = 0;
            r_cross_X_z(smtheta) = 0;
            X_x(bgtheta) = 0;
            X_y(bgtheta) = 0;
            X_z(bgtheta) = 0;
            r_cross_X_x(bgtheta) = 0;
            r_cross_X_y(bgtheta) = 0;
            r_cross_X_z(bgtheta) = 0;
        end
    else % RL: at the moment only for theta ~ 0, need to do theta ~ pi
        csc_Y_l_m = cscfieldtheta.*Y_l_m;
        csc_Y_lminus1_m = cscfieldtheta.*Y_lminus1_m;
        if m == -1
            csc_Y_l_m(smtheta) = csc_Y_l_m_series_expans(fieldtheta(smtheta),fieldphi(smtheta),l,m);
            csc_Y_lminus1_m(smtheta) = csc_Y_l_m_series_expans(fieldtheta(smtheta),fieldphi(smtheta),l-1,m);
        elseif m == 1
            csc_Y_l_m(smtheta) = -conj(csc_Y_l_m_series_expans(fieldtheta(smtheta),fieldphi(smtheta),l,-m));
            csc_Y_lminus1_m(smtheta) = -conj(csc_Y_l_m_series_expans(fieldtheta(smtheta),fieldphi(smtheta),l-1,-m));
        else
            csc_Y_l_m(smtheta) = 0;
            csc_Y_lminus1_m(smtheta) = 0;
        end

            X_x(smtheta) = (1/lnorm).*((-m*cosfieldphi(smtheta)+1i*l*sinfieldphi(smtheta)).*cosfieldtheta(smtheta).*csc_Y_l_m(smtheta) ...
                - 1i*lmul.*sinfieldphi(smtheta).*csc_Y_lminus1_m(smtheta));
            X_y(smtheta) = (1/lnorm).*((-m*sinfieldphi(smtheta)-1i*l*cosfieldphi(smtheta)).*cosfieldtheta(smtheta).*csc_Y_l_m(smtheta) ...
                + 1i*lmul.*cosfieldphi(smtheta).*csc_Y_lminus1_m(smtheta));
            X_z(smtheta) = (1/lnorm).*m.*Y_l_m(smtheta);
            
            r_cross_X_x(smtheta) = (1/lnorm).*((m*sinfieldphi(smtheta)+1i*l*cos2fieldtheta(smtheta).*cosfieldphi(smtheta)).*csc_Y_l_m(smtheta) ...
                - 1i*lmul*cosfieldtheta(smtheta).*cosfieldphi(smtheta).*csc_Y_lminus1_m(smtheta));
            r_cross_X_y(smtheta) = (1/lnorm).*((-m*cosfieldphi(smtheta)+1i*l*cos2fieldtheta(smtheta).*sinfieldphi(smtheta)).*csc_Y_l_m(smtheta) ...
                - 1i*lmul*cosfieldtheta(smtheta).*sinfieldphi(smtheta).*csc_Y_lminus1_m(smtheta));
            r_cross_X_z(smtheta) = (-1i/lnorm).*(l*cosfieldtheta(smtheta).*csc_Y_l_m(smtheta) ...
                - lmul.*Y_lminus1_m(smtheta));
            
            X_x(bgtheta) = (1/lnorm).*((-m*cosfieldphi(bgtheta)+1i*l*sinfieldphi(bgtheta)).*cosfieldtheta(bgtheta).*csc_Y_l_m(bgtheta) ...
                - 1i*lmul.*sinfieldphi(bgtheta).*csc_Y_lminus1_m(bgtheta));
            X_y(bgtheta) = (1/lnorm).*((-m*sinfieldphi(bgtheta)-1i*l*cosfieldphi(bgtheta)).*cosfieldtheta(bgtheta).*csc_Y_l_m(bgtheta) ...
                + 1i*lmul.*cosfieldphi(bgtheta).*csc_Y_lminus1_m(bgtheta));
            X_z(bgtheta) = (1/lnorm).*m.*Y_l_m(bgtheta);
            
            r_cross_X_x(bgtheta) = (1/lnorm).*((m*sinfieldphi(bgtheta)+1i*l*cos2fieldtheta(bgtheta).*cosfieldphi(bgtheta)).*csc_Y_l_m(bgtheta) ...
                - 1i*lmul*cosfieldtheta(bgtheta).*cosfieldphi(bgtheta).*csc_Y_lminus1_m(bgtheta));
            r_cross_X_y(bgtheta) = (1/lnorm).*((-m*cosfieldphi(bgtheta)+1i*l*cos2fieldtheta(bgtheta).*sinfieldphi(bgtheta)).*csc_Y_l_m(bgtheta) ...
                - 1i*lmul*cosfieldtheta(bgtheta).*sinfieldphi(bgtheta).*csc_Y_lminus1_m(bgtheta));
            r_cross_X_z(bgtheta) = (-1i/lnorm).*(l*cosfieldtheta(bgtheta).*csc_Y_l_m(bgtheta) ...
                - lmul.*Y_lminus1_m(bgtheta));
    end
end

M_x = J_l_kr.*X_x;
M_y = J_l_kr.*X_y;
M_z = J_l_kr.*X_z;
M3_x = H_l_1_kr.*X_x;
M3_y = H_l_1_kr.*X_y;
M3_z = H_l_1_kr.*X_z;

N_x = (1./krset).*J_l_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_x;
N_y = (1./krset).*J_l_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_y;
N_z = (1./krset).*J_l_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_z;

N3_x = (1./krset).*H_l_1_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_x;
N3_y = (1./krset).*H_l_1_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_y;
N3_z = (1./krset).*H_l_1_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_z;

smrho = find(fieldrho < eps);
M_x(smrho) = 0;
M_y(smrho) = 0;
M_z(smrho) = 0;
M3_x(smrho) = 0;
M3_y(smrho) = 0;
M3_z(smrho) = 0;
N3_x(smrho) = 0; % It doesn't matter since these are not used at the origin
N3_y(smrho) = 0;
N3_z(smrho) = 0;
if l == 1,
    switch m
        case -1
            N_x(smrho)=1i*sqrt(1/(6*pi))*(1/sqrt(2));
            N_y(smrho)=1i*sqrt(1/(6*pi))*(-1i/sqrt(2));
            N_z(smrho)=0;
        case 0
            N_x(smrho)=0;
            N_y(smrho)=0;
            N_z(smrho)=1i*sqrt(1/(6*pi));
        case 1
            N_x(smrho)=1i*sqrt(1/(6*pi))*(-1/sqrt(2));
            N_y(smrho)=1i*sqrt(1/(6*pi))*(-1i/sqrt(2));
            N_z(smrho)=0;
    end
else
    N_x(smrho)=0;
    N_y(smrho)=0;
    N_z(smrho)=0;
end

% if fieldrho < eps
%     M_x = 0;
%     M_y = 0;
%     M_z = 0;
%     M3_x = 0;
%     M3_y = 0;
%     M3_z = 0;
%     N3_x=0; % It doesn't matter since these are not used at the origin
%     N3_y=0;
%     N3_z=0;
%     if l == 1,
%         switch m
%             case -1
%                 N_x=1i*sqrt(1/(6*pi))*(1/sqrt(2));
%                 N_y=1i*sqrt(1/(6*pi))*(-1i/sqrt(2));
%                 N_z=0;
%             case 0
%                 N_x=0;
%                 N_y=0;
%                 N_z=1i*sqrt(1/(6*pi));
%             case 1
%                 N_x=1i*sqrt(1/(6*pi))*(-1/sqrt(2));
%                 N_y=1i*sqrt(1/(6*pi))*(-1i/sqrt(2));
%                 N_z=0;
%         end
%     else
%         N_x=0;
%         N_y=0;
%         N_z=0;
%         N3_x=0;
%         N3_y=0;
%         N3_z=0;
%     end
% else
%     M_x = J_l_kr.*X_x;
%     M_y = J_l_kr.*X_y;
%     M_z = J_l_kr.*X_z;
%     M3_x = H_l_1_kr.*X_x;
%     M3_y = H_l_1_kr.*X_y;
%     M3_z = H_l_1_kr.*X_z;
%     
%     N_x = (1./krset).*J_l_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_x;
%     N_y = (1./krset).*J_l_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_y;
%     N_z = (1./krset).*J_l_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*J_l_kr.*Y_l_m.*rhat_z;
%     
%     N3_x = (1./krset).*H_l_1_kr_prime.*r_cross_X_x + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_x;
%     N3_y = (1./krset).*H_l_1_kr_prime.*r_cross_X_y + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_y;
%     N3_z = (1./krset).*H_l_1_kr_prime.*r_cross_X_z + 1i*(lnorm./krset).*H_l_1_kr.*Y_l_m.*rhat_z;
% end

T_hat_field = T_matrix_set{counter_fields};

E_x(:,mask_r4) = ele_scaling*(T_hat_field*[ (D_m_41(mask_r4).*M_x(mask_r4)).'; (D_n_41(mask_r4).*N_x(mask_r4)).' ]);
E_y(:,mask_r4) = ele_scaling*(T_hat_field*[ (D_m_41(mask_r4).*M_y(mask_r4)).'; (D_n_41(mask_r4).*N_y(mask_r4)).' ]);
E_z(:,mask_r4) = ele_scaling*(T_hat_field*[ (D_m_41(mask_r4).*M_z(mask_r4)).'; (D_n_41(mask_r4).*N_z(mask_r4)).' ]);
B_x(:,mask_r4) = mag_scaling_r4*(T_hat_field*[ (D_m_41(mask_r4).*N_x(mask_r4)).'; (D_n_41(mask_r4).*M_x(mask_r4)).' ]);
B_y(:,mask_r4) = mag_scaling_r4*(T_hat_field*[ (D_m_41(mask_r4).*N_y(mask_r4)).'; (D_n_41(mask_r4).*M_y(mask_r4)).' ]);
B_z(:,mask_r4) = mag_scaling_r4*(T_hat_field*[ (D_m_41(mask_r4).*N_z(mask_r4)).'; (D_n_41(mask_r4).*M_z(mask_r4)).' ]);

E_x(:,mask_r3) = ele_scaling*(T_hat_field*[ (B_m_31(mask_r3).*M3_x(mask_r3) + D_m_31(mask_r3).*M_x(mask_r3)).'; (B_n_31(mask_r3).*N3_x(mask_r3) + D_n_31(mask_r3).*N_x(mask_r3)).' ]);
E_y(:,mask_r3) = ele_scaling*(T_hat_field*[ (B_m_31(mask_r3).*M3_y(mask_r3) + D_m_31(mask_r3).*M_y(mask_r3)).'; (B_n_31(mask_r3).*N3_y(mask_r3) + D_n_31(mask_r3).*N_y(mask_r3)).' ]);
E_z(:,mask_r3) = ele_scaling*(T_hat_field*[ (B_m_31(mask_r3).*M3_z(mask_r3) + D_m_31(mask_r3).*M_z(mask_r3)).'; (B_n_31(mask_r3).*N3_z(mask_r3) + D_n_31(mask_r3).*N_z(mask_r3)).' ]);
B_x(:,mask_r3) = mag_scaling_r3*(T_hat_field*[ (B_m_31(mask_r3).*N3_x(mask_r3) + D_m_31(mask_r3).*N_x(mask_r3)).'; (B_n_31(mask_r3).*M3_x(mask_r3) + D_n_31(mask_r3).*M_x(mask_r3)).' ]);
B_y(:,mask_r3) = mag_scaling_r3*(T_hat_field*[ (B_m_31(mask_r3).*N3_y(mask_r3) + D_m_31(mask_r3).*N_y(mask_r3)).'; (B_n_31(mask_r3).*M3_y(mask_r3) + D_n_31(mask_r3).*M_y(mask_r3)).' ]);
B_z(:,mask_r3) = mag_scaling_r3*(T_hat_field*[ (B_m_31(mask_r3).*N3_z(mask_r3) + D_m_31(mask_r3).*N_z(mask_r3)).'; (B_n_31(mask_r3).*M3_z(mask_r3) + D_n_31(mask_r3).*M_z(mask_r3)).' ]);

E_x(:,mask_r2) = ele_scaling*(T_hat_field*[ (B_m_21(mask_r2).*M3_x(mask_r2) + D_m_21(mask_r2).*M_x(mask_r2)).'; (B_n_21(mask_r2).*N3_x(mask_r2) + D_n_21(mask_r2).*N_x(mask_r2)).' ]);
E_y(:,mask_r2) = ele_scaling*(T_hat_field*[ (B_m_21(mask_r2).*M3_y(mask_r2) + D_m_21(mask_r2).*M_y(mask_r2)).'; (B_n_21(mask_r2).*N3_y(mask_r2) + D_n_21(mask_r2).*N_y(mask_r2)).' ]);
E_z(:,mask_r2) = ele_scaling*(T_hat_field*[ (B_m_21(mask_r2).*M3_z(mask_r2) + D_m_21(mask_r2).*M_z(mask_r2)).'; (B_n_21(mask_r2).*N3_z(mask_r2) + D_n_21(mask_r2).*N_z(mask_r2)).' ]);
B_x(:,mask_r2) = mag_scaling_r2*(T_hat_field*[ (B_m_21(mask_r2).*N3_x(mask_r2) + D_m_21(mask_r2).*N_x(mask_r2)).'; (B_n_21(mask_r2).*M3_x(mask_r2) + D_n_21(mask_r2).*M_x(mask_r2)).' ]);
B_y(:,mask_r2) = mag_scaling_r2*(T_hat_field*[ (B_m_21(mask_r2).*N3_y(mask_r2) + D_m_21(mask_r2).*N_y(mask_r2)).'; (B_n_21(mask_r2).*M3_y(mask_r2) + D_n_21(mask_r2).*M_y(mask_r2)).' ]);
B_z(:,mask_r2) = mag_scaling_r2*(T_hat_field*[ (B_m_21(mask_r2).*N3_z(mask_r2) + D_m_21(mask_r2).*N_z(mask_r2)).'; (B_n_21(mask_r2).*M3_z(mask_r2) + D_n_21(mask_r2).*M_z(mask_r2)).' ]);

%         E_x(:,mask_r1) = TBD;
%         E_y(:,mask_r1) = TBD;
%         E_z(:,mask_r1) = TBD;
%         B_x(:,mask_r1) = TBD;
%         B_y(:,mask_r1) = TBD;
%         B_z(:,mask_r1) = TBD;


% form field contribution for the current mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(whichcurrents)==1,
    fieldweights = currentweights(counter_fields);
else
    fieldweights = currentweights(counter_fields,:);
end
Efpt(:,:,:,1) = reshape(fieldweights*reshape(E_x,[2 numel(fieldtheta)]),size(fieldtheta));
Efpt(:,:,:,2) = reshape(fieldweights*reshape(E_y,[2 numel(fieldtheta)]),size(fieldtheta));
Efpt(:,:,:,3) = reshape(fieldweights*reshape(E_z,[2 numel(fieldtheta)]),size(fieldtheta));
Bfpt(:,:,:,1) = reshape(fieldweights*reshape(B_x,[2 numel(fieldtheta)]),size(fieldtheta));
Bfpt(:,:,:,2) = reshape(fieldweights*reshape(B_y,[2 numel(fieldtheta)]),size(fieldtheta));
Bfpt(:,:,:,3) = reshape(fieldweights*reshape(B_z,[2 numel(fieldtheta)]),size(fieldtheta));

% add weighted fields to running pattern buffer %%%%%%%%%%%%%%%%%%%%%%%%%%%
Efieldpattern(:,:,:,:,find(test_1 & test_2)) = Efieldpattern(:,:,:,:,find(test_1 & test_2)) + Efpt;
Bfieldpattern(:,:,:,:,find(test_1 & test_2)) = Bfieldpattern(:,:,:,:,find(test_1 & test_2)) + Bfpt;


% 
% 
% 
% 
% 
% 
% % coil field pattern %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if s_opts.compute_coil_current_pattern
%     if (~include_eldip) && include_cylwind
%         W_coil = reshape(W_cylwind,[nm nn ncoil]);
%         if p_opts.unit_current_weighting
%             netcurrentweights_coil = squeeze(W_coil(find(m==whichm),find(n==whichn),:)).'*ones(ncoil,1);
%         else
%             netcurrentweights_coil = squeeze(W_coil(find(m==whichm),find(n==whichn),:)).'*currentweights_coil(:);
%         end
%         Efpt_coil(:,:,:,1) = reshape(netcurrentweights_coil*E_x(1,:),size(fieldtheta));
%         Efpt_coil(:,:,:,2) = reshape(netcurrentweights_coil*E_y(1,:),size(fieldtheta));
%         Efpt_coil(:,:,:,3) = reshape(netcurrentweights_coil*E_z(1,:),size(fieldtheta));
%         Bfpt_coil(:,:,:,1) = reshape(netcurrentweights_coil*B_x(1,:),size(fieldtheta));
%         Bfpt_coil(:,:,:,2) = reshape(netcurrentweights_coil*B_y(1,:),size(fieldtheta));
%         Bfpt_coil(:,:,:,3) = reshape(netcurrentweights_coil*B_z(1,:),size(fieldtheta));
%         Efieldpattern_coil(:,:,:,:,find(test_1 & test_2)) = Efieldpattern_coil(:,:,:,:,find(test_1 & test_2)) + Efpt_coil;
%         Bfieldpattern_coil(:,:,:,:,find(test_1 & test_2)) = Bfieldpattern_coil(:,:,:,:,find(test_1 & test_2)) + Bfpt_coil;
%     elseif include_eldip && (~include_cylwind)
%         W_coil_mag = reshape(W_eldip(1:2:numbasis,:),[nm nn ndipole]);
%         W_coil_ele = reshape(W_eldip(2:2:numbasis,:),[nm nn ndipole]);
%         if p_opts.unit_current_weighting
%             netcurrentweights_coil_M = squeeze(W_coil_mag(find(m==whichm),find(n==whichn),:)).'*ones(ndipole,1);
%             netcurrentweights_coil_E = squeeze(W_coil_ele(find(m==whichm),find(n==whichn),:)).'*ones(ndipole,1);
%         else
%             netcurrentweights_coil_M = squeeze(W_coil_mag(find(m==whichm),find(n==whichn),1,:)).'*currentweights_coil(:);
%             netcurrentweights_coil_E = squeeze(W_coil_ele(find(m==whichm),find(n==whichn),2,:)).'*currentweights_coil(:);
%         end
%         netcurrentweights_coil = [netcurrentweights_coil_M netcurrentweights_coil_E];
%         Efpt_coil(:,:,:,1) = reshape(netcurrentweights_coil*E_x,size(fieldtheta));
%         Efpt_coil(:,:,:,2) = reshape(netcurrentweights_coil*E_y,size(fieldtheta));
%         Efpt_coil(:,:,:,3) = reshape(netcurrentweights_coil*E_z,size(fieldtheta));
%         Bfpt_coil(:,:,:,1) = reshape(netcurrentweights_coil*B_x,size(fieldtheta));
%         Bfpt_coil(:,:,:,2) = reshape(netcurrentweights_coil*B_y,size(fieldtheta));
%         Bfpt_coil(:,:,:,3) = reshape(netcurrentweights_coil*B_z,size(fieldtheta));
%         Efieldpattern_coil(:,:,:,:,find(test_1 & test_2)) = Efieldpattern_coil(:,:,:,:,find(test_1 & test_2)) + Efpt_coil;
%         Bfieldpattern_coil(:,:,:,:,find(test_1 & test_2)) = Bfieldpattern_coil(:,:,:,:,find(test_1 & test_2)) + Bfpt_coil;        
%         if p_opts.dipole_curr_decomposition
%             Efpt_coil_mag(:,:,:,1) = reshape(netcurrentweights_coil(1)*E_x(1,:),size(fieldtheta));
%             Efpt_coil_mag(:,:,:,2) = reshape(netcurrentweights_coil(1)*E_y(1,:),size(fieldtheta));
%             Efpt_coil_mag(:,:,:,3) = reshape(netcurrentweights_coil(1)*E_z(1,:),size(fieldtheta));
%             Bfpt_coil_mag(:,:,:,1) = reshape(netcurrentweights_coil(1)*B_x(1,:),size(fieldtheta));
%             Bfpt_coil_mag(:,:,:,2) = reshape(netcurrentweights_coil(1)*B_y(1,:),size(fieldtheta));
%             Bfpt_coil_mag(:,:,:,3) = reshape(netcurrentweights_coil(1)*B_z(1,:),size(fieldtheta));
%             Efieldpattern_dipole_mag(:,:,:,:,find(test_1 & test_2)) = Efieldpattern_dipole_mag(:,:,:,:,find(test_1 & test_2)) + Efpt_coil_mag;
%             Bfieldpattern_dipole_mag(:,:,:,:,find(test_1 & test_2)) = Bfieldpattern_dipole_mag(:,:,:,:,find(test_1 & test_2)) + Bfpt_coil_mag;
%             Efpt_coil_ele(:,:,:,1) = reshape(netcurrentweights_coil(2)*E_x(2,:),size(fieldtheta));
%             Efpt_coil_ele(:,:,:,2) = reshape(netcurrentweights_coil(2)*E_y(2,:),size(fieldtheta));
%             Efpt_coil_ele(:,:,:,3) = reshape(netcurrentweights_coil(2)*E_z(2,:),size(fieldtheta));
%             Bfpt_coil_ele(:,:,:,1) = reshape(netcurrentweights_coil(2)*B_x(2,:),size(fieldtheta));
%             Bfpt_coil_ele(:,:,:,2) = reshape(netcurrentweights_coil(2)*B_y(2,:),size(fieldtheta));
%             Bfpt_coil_ele(:,:,:,3) = reshape(netcurrentweights_coil(2)*B_z(2,:),size(fieldtheta));
%             Efieldpattern_dipole_ele(:,:,:,:,find(test_1 & test_2)) = Efieldpattern_dipole_ele(:,:,:,:,find(test_1 & test_2)) + Efpt_coil_ele;
%             Bfieldpattern_dipole_ele(:,:,:,:,find(test_1 & test_2)) = Bfieldpattern_dipole_ele(:,:,:,:,find(test_1 & test_2)) + Bfpt_coil_ele;
%         end
%     else % both cylindrical window coils and electric dipoles
%         W_cylwind_mag = reshape(W_cylwind,[nm nn ncoil]);
%         W_dipole_mag = reshape(W_eldip(1:2:numbasis,:),[nm nn ndipole]);
%         W_dipole_ele = reshape(W_eldip(2:2:numbasis,:),[nm nn ndipole]);
%         
%         currentweights_coil_cylwind = currentweights_coil(1:ncoil,:);
%         if p_opts.unit_current_weighting
%             netcurrentweights_cylwind = squeeze(W_cylwind_mag(find(m==whichm),find(n==whichn),:)).'*ones(ncoil,1);
%         else
%             netcurrentweights_cylwind = squeeze(W_cylwind_mag(find(m==whichm),find(n==whichn),:)).'*currentweights_coil_cylwind(:);
%         end
%         
%         currentweights_coil_dipole = currentweights_coil(ncoil+1:ncoil+ndipole,:);
%         if p_opts.unit_current_weighting
%             netcurrentweights_dipole_M = squeeze(W_dipole_mag(find(m==whichm),find(n==whichn),:)).'*ones(ndipole,1);
%             netcurrentweights_dipole_E = squeeze(W_dipole_ele(find(m==whichm),find(n==whichn),:)).'*ones(ndipole,1);
%         else
%             netcurrentweights_dipole_M = squeeze(W_dipole_mag(find(m==whichm),find(n==whichn),1,:)).'*currentweights_coil_dipole(:);
%             netcurrentweights_dipole_E = squeeze(W_dipole_ele(find(m==whichm),find(n==whichn),2,:)).'*currentweights_coil_dipole(:);
%         end
%         netcurrentweights_dipole = [netcurrentweights_dipole_M netcurrentweights_dipole_E];
%         
%         Efpt_coil(:,:,:,1) = reshape(netcurrentweights_cylwind*E_x(1,:)+netcurrentweights_dipole*E_x,size(fieldtheta));
%         Efpt_coil(:,:,:,2) = reshape(netcurrentweights_cylwind*E_y(1,:)+netcurrentweights_dipole*E_y,size(fieldtheta));
%         Efpt_coil(:,:,:,3) = reshape(netcurrentweights_cylwind*E_z(1,:)+netcurrentweights_dipole*E_z,size(fieldtheta));
%         Bfpt_coil(:,:,:,1) = reshape(netcurrentweights_cylwind*B_x(1,:)+netcurrentweights_dipole*B_x,size(fieldtheta));
%         Bfpt_coil(:,:,:,2) = reshape(netcurrentweights_cylwind*B_y(1,:)+netcurrentweights_dipole*B_y,size(fieldtheta));
%         Bfpt_coil(:,:,:,3) = reshape(netcurrentweights_cylwind*B_z(1,:)+netcurrentweights_dipole*B_z,size(fieldtheta));
%         Efieldpattern_coil(:,:,:,:,find(test_1 & test_2)) = Efieldpattern_coil(:,:,:,:,find(test_1 & test_2)) + Efpt_coil;
%         Bfieldpattern_coil(:,:,:,:,find(test_1 & test_2)) = Bfieldpattern_coil(:,:,:,:,find(test_1 & test_2)) + Bfpt_coil;        
%         
%     end
% end