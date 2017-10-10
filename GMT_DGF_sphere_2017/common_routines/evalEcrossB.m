function [ExB] = evalEcrossB(r,theta,phi,whichcurrents,sphereradius,k_0,k_0_rad_outer,k_s,l,m,ele_scaling,mag_scaling)

% xset=r*cos(phi)*sin(theta);
% yset= r*sin(phi)*sin(theta);
% zset=r*cos(theta);

costhetaset = cos(theta);
sinthetaset = sin(theta);                % sin[theta]
sinphiset = sin(phi);                    % sin[phi]
cosphiset = cos(phi);                    % cos[phi]

cotthetaset = cot(theta);                % cot[theta] ** Inf at theta=0
cscthetaset = csc(theta);                % csc[theta] ** Inf at theta=0
cos2thetaset = costhetaset.*costhetaset; % cos[theta]^2


% unit vectors conversion
rhat_x = sinthetaset.*cosphiset;
rhat_y = sinthetaset.*sinphiset;
rhat_z = costhetaset;

% thetahat_x = costhetaset.*cosphiset;
% thetahat_y = costhetaset.*sinphiset;
% thetahat_z = -sinthetaset;
% 
% phihat_x = -sinphiset;
% phihat_y = cosphiset;

k_0r = k_0*r;
besselnorm = sqrt(pi./(2.*k_0r));
lnorm = sqrt(l*(l+1));
lmul = sqrt((2*l + 1)*(l^2 - m^2)/(2*l - 1));
besselnorm_rad_outer = sqrt(pi/(2*k_0_rad_outer));

% ----------------------------------------- %
H_l_1_k_0r = besselnorm*besselh(l+0.5,1,k_0r);
H_l_1_k_0r_prime = -k_0r*besselnorm*besselh(l+1+0.5,1,k_0r) + (l + 1)*H_l_1_k_0r; % this is actually d(r*hl)/dr

% J_l_k_0r = besselnorm.*besselj(l+0.5,k_0r); % DEFINING BESSEL FOR K0 at rset ( poynting radius)
% J_l_k_0r(r<eps) = (l==0);
% J_lplus1_k_0r = besselnorm.*besselj(l+1+0.5,k_0r);
% J_lplus1_k_0r(r<eps) = 0;
% J_l_k_0r_prime = -k_0r.*J_lplus1_k_0r + (l + 1)*J_l_k_0r;


H_l_1_k_0_rad_outer = besselnorm_rad_outer*besselh(l+0.5,1,k_0_rad_outer);
H_l_1_k_0_rad_outer_prime = -k_0_rad_outer*besselnorm_rad_outer*besselh(l+1+0.5,1,k_0_rad_outer) + (l + 1)*H_l_1_k_0_rad_outer;

J_l_k_0_outer_radius = besselnorm_rad_outer.*besselj(l+0.5,k_0_rad_outer); % DEFINING BESSEL FOR K0 at r' ( outer_radius)
J_l_k_0_outer_radius(r<eps) = (l==0);
J_lplus1_k_0_outer_radius = besselnorm_rad_outer.*besselj(l+1+0.5,k_0_rad_outer);
J_lplus1_k_0_outer_radius(r<eps) = 0;
J_l_k_0_outer_radius_prime = -k_0_rad_outer.*J_lplus1_k_0_outer_radius + (l + 1)*J_l_k_0_outer_radius;


[a_l, b_l, c_l, d_l] = compute_T_coefficients_sphere(l,k_0,k_s,sphereradius);


T = -1i*lnorm*([J_l_k_0_outer_radius + H_l_1_k_0_rad_outer*a_l                    0
    0                 (1/k_0_rad_outer)*(J_l_k_0_outer_radius_prime + H_l_1_k_0_rad_outer_prime*b_l)].');% changed from rset to outer_radius  % for j, do we need to convert to cartesian? 2/7/2012



legendrenorm = sqrt((2*l + 1)/(4*pi));
legendrenorm_minus1 = sqrt((2*l - 1)/(4*pi));
legendrefunctions = legendre(l,costhetaset,'sch'); % # row is m = 0,...l ; # col is R
if l>0,
    legendrefunctions_lminus1 = legendre(l-1,costhetaset,'sch');
    legendrefunctions_lminus1 = [legendrefunctions_lminus1; zeros(size(costhetaset))];
end
if m > 0
    Y_l_m = ((-1)^m)*legendrenorm*(1/sqrt(2))*legendrefunctions((m+1),:).*exp(1i*m*phi);
elseif m == 0
    Y_l_m = legendrenorm*legendrefunctions(1,:);
else % using Y_(l,-m) = (-1)^m*conj[Y_(l,m)]
    Y_l_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm*(1/sqrt(2))*legendrefunctions((abs(m)+1),:).*exp(1i*abs(m)*phi));
end

if m > 0
    Y_lminus1_m = ((-1)^m)*legendrenorm_minus1*(1/sqrt(2))*legendrefunctions_lminus1(m+1,:).*exp(1i*m*phi);
elseif m == 0
    Y_lminus1_m = legendrenorm_minus1*legendrefunctions_lminus1(1,:);
else
    Y_lminus1_m = ((-1)^abs(m))*conj(((-1)^abs(m))*legendrenorm_minus1*(1/sqrt(2))*legendrefunctions_lminus1((abs(m)+1),:).*exp(i*abs(m)*phi));
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
    - 1i*lmul.*Y_lminus1_m);


T_hat = T;
if whichcurrents == 1,
    T_hat = [1 0]*T;
end
if whichcurrents == 2,
    T_hat = [0 1]*T;
end


% M_x = J_l_k_0r.*X_x;  % rset is at poynting surface r
% M_y = J_l_k_0r.*X_y;
% M_z = J_l_k_0r.*X_z;

% CALCULATING M+--------------
M_x_H =  H_l_1_k_0r.*X_x; % CHANGED TO HANKEL FUNCTION
M_y_H =  H_l_1_k_0r.*X_y; % THIS IS EQUAL TO THE M+ IN THESIS
M_z_H =  H_l_1_k_0r.*X_z; % ALL MS AND NS ARE M~ AND N~


%********  Verified calculation **********
N_x_H = (1/k_0r).* H_l_1_k_0r_prime.*r_cross_X_x + 1i*(lnorm./k_0r).* H_l_1_k_0r.*Y_l_m.*rhat_x;
N_y_H = (1/k_0r).* H_l_1_k_0r_prime.*r_cross_X_y + 1i*(lnorm./k_0r).* H_l_1_k_0r.*Y_l_m.*rhat_y;
N_z_H = (1/k_0r).* H_l_1_k_0r_prime.*r_cross_X_z + 1i*(lnorm./k_0r).* H_l_1_k_0r.*Y_l_m.*rhat_z;


%********  STORE ALL FIELDS **********
E_x = ele_scaling*(T_hat*[M_x_H; N_x_H]);
E_y = ele_scaling*(T_hat*[M_y_H; N_y_H]);
E_z = ele_scaling*(T_hat*[M_z_H; N_z_H]);

B_x = mag_scaling*(T_hat*[N_x_H; M_x_H]);
B_y = mag_scaling*(T_hat*[N_y_H; M_y_H]);
B_z = mag_scaling*(T_hat*[N_z_H; M_z_H]);

E = [E_x E_y E_z];
B = [B_x B_y B_z];
ExB(1,:) = cross(E(1,:), conj(B(1,:)));
ExB(2,:) = cross(E(2,:), conj(B(2,:)));
end % end function ExB




