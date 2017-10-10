function    [R_H_P_f, R_H_F_f, R_V_P_f, R_V_F_f, T_H_P_f, T_H_F_f, T_V_P_f, T_V_F_f] = ...
    compute_reflection_and_transmission_coef(l,k_f,k_fplus1,rad_f)

%
% taken from: Li et al. 
%             Electromagnetic Dyadic Green's Function in Spherically Multilayered Media
%             IEEE Transactions on Microwave Theory and Techniques, vol. 42(12), 1994
%
% input:
%   l: mode index
%   k_f: wavenumber in f-th layer (m^-1)
%   k_fplus1: wavenumber in (f + 1)-th layer (m^-1)
%   rad_f: radius separating layer f from (f + 1) (m)
%
% output:
%   R_H_P_f, R_H_F_f, R_V_P_f, R_V_F_f, T_H_P_f, T_H_F_f, T_V_P_f, T_V_F_f: values of the appropriate coefficients
%
% Riccardo Lattanzi, 9/3/13
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 4*pi*1e-7;         % permeability of free space [Wb][A^-1][m^-1]

krad_f_f = k_f*rad_f;
krad_fplus1_f = k_fplus1*rad_f;

besselnorm_f_f = sqrt(pi/(2*krad_f_f));
besselnorm_fplus1_f = sqrt(pi/(2*krad_fplus1_f));

J_l_krad_f_f = besselnorm_f_f*besselj(l+0.5,krad_f_f);
J_l_krad_fplus1_f = besselnorm_fplus1_f*besselj(l+0.5,krad_fplus1_f);

J_l_krad_f_f_prime = -besselnorm_f_f*besselj(l+1+0.5,krad_f_f) + ((l + 1)/krad_f_f)*J_l_krad_f_f;
J_l_krad_fplus1_f_prime = -besselnorm_fplus1_f*besselj(l+1+0.5,krad_fplus1_f) + ((l + 1)/krad_fplus1_f)*J_l_krad_fplus1_f;


H_l_1_krad_f_f = besselnorm_f_f*besselh(l+0.5,1,krad_f_f);
H_l_1_krad_fplus1_f = besselnorm_fplus1_f*besselh(l+0.5,1,krad_fplus1_f);

% H_l_1_krad_f_f_prime = -besselh(l+1+0.5,1,krad_f_f); + ((l + 1)/krad_f_f)*H_l_1_krad_f_f;
% H_l_1_krad_fplus1_f_prime = -besselh(l+1+0.5,1,krad_fplus1_f); + ((l + 1)/krad_fplus1_f)*H_l_1_krad_fplus1_f;

% RL Dec 27, 2013 removed semicolon
H_l_1_krad_f_f_prime = -besselnorm_f_f*besselh(l+1+0.5,1,krad_f_f) + ((l + 1)/krad_f_f)*H_l_1_krad_f_f;
H_l_1_krad_fplus1_f_prime = -besselnorm_fplus1_f*besselh(l+1+0.5,1,krad_fplus1_f) + ((l + 1)/krad_fplus1_f)*H_l_1_krad_fplus1_f;


R_H_P_f = (mu*k_fplus1*H_l_1_krad_fplus1_f_prime*H_l_1_krad_f_f - mu*k_f*H_l_1_krad_f_f_prime*H_l_1_krad_fplus1_f)/...
              (mu*k_fplus1*J_l_krad_f_f*H_l_1_krad_fplus1_f_prime - mu*k_f*J_l_krad_f_f_prime*H_l_1_krad_fplus1_f);
          
R_H_F_f = (mu*k_fplus1*J_l_krad_fplus1_f_prime*J_l_krad_f_f - mu*k_f*J_l_krad_f_f_prime*J_l_krad_fplus1_f)/...
              (mu*k_fplus1*J_l_krad_fplus1_f_prime*H_l_1_krad_f_f - mu*k_f*J_l_krad_fplus1_f*H_l_1_krad_f_f_prime);
          
R_V_P_f = (mu*k_fplus1*H_l_1_krad_fplus1_f*H_l_1_krad_f_f_prime - mu*k_f*H_l_1_krad_f_f*H_l_1_krad_fplus1_f_prime)/...
              (mu*k_fplus1*J_l_krad_f_f_prime*H_l_1_krad_fplus1_f - mu*k_f*J_l_krad_f_f*H_l_1_krad_fplus1_f_prime);
          
R_V_F_f = (mu*k_fplus1*J_l_krad_fplus1_f*J_l_krad_f_f_prime - mu*k_f*J_l_krad_f_f*J_l_krad_fplus1_f_prime)/...
              (mu*k_fplus1*J_l_krad_fplus1_f*H_l_1_krad_f_f_prime - mu*k_f*J_l_krad_fplus1_f_prime*H_l_1_krad_f_f);
          
          
T_H_P_f = mu*k_fplus1*(J_l_krad_fplus1_f*H_l_1_krad_fplus1_f_prime - J_l_krad_fplus1_f_prime*H_l_1_krad_fplus1_f)/...
              (mu*k_fplus1*J_l_krad_f_f*H_l_1_krad_fplus1_f_prime - mu*k_f*J_l_krad_f_f_prime*H_l_1_krad_fplus1_f);

T_H_F_f = mu*k_fplus1*(J_l_krad_fplus1_f_prime*H_l_1_krad_fplus1_f - J_l_krad_fplus1_f*H_l_1_krad_fplus1_f_prime)/...
              (mu*k_fplus1*J_l_krad_fplus1_f_prime*H_l_1_krad_f_f - mu*k_f*J_l_krad_fplus1_f*H_l_1_krad_f_f_prime);
          
T_V_P_f = mu*k_fplus1*(J_l_krad_fplus1_f_prime*H_l_1_krad_fplus1_f - J_l_krad_fplus1_f*H_l_1_krad_fplus1_f_prime)/...
              (mu*k_fplus1*J_l_krad_f_f_prime*H_l_1_krad_fplus1_f - mu*k_f*J_l_krad_f_f*H_l_1_krad_fplus1_f_prime);
          
T_V_F_f = mu*k_fplus1*(J_l_krad_fplus1_f*H_l_1_krad_fplus1_f_prime - J_l_krad_fplus1_f_prime*H_l_1_krad_fplus1_f)/...
              (mu*k_fplus1*J_l_krad_fplus1_f*H_l_1_krad_f_f_prime - mu*k_f*J_l_krad_fplus1_f_prime*H_l_1_krad_f_f);
          
          
          
          



