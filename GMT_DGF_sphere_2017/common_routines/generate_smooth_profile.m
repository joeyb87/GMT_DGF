% 2D Gaussian profile

function [gaussian_profile] = generate_smooth_profile(sph_rad, X, Y, mask, plot_flag, prof_type)

A = 1; % central peak
x0 = 0; y0 = 0;
sigma_x = 0.1;
sigma_y = 0.1;

Z = A*exp(-( ((X-x0)./sigma_x).^2 + ((Y-y0)./sigma_y).^2 ));

switch prof_type
    case 1
        gaussian_profile = Z;
        gaussian_profile(~mask) = 0;
    case 2
        Z_temp = -Z;
        gaussian_profile = Z_temp -min(Z_temp(:));
        gaussian_profile(~mask) = 0;
end

if plot_flag
    figure;
    set(gcf,'name','Target Excitation Profile in 3D');
    surf(X,Y,gaussian_profile);
    title(['Target Excitation Profile'])
    iptsetpref('ImshowInitialMagnification','fit');
    figure;
    set(gcf,'name','Target Excitation Profile');
    imshow(gaussian_profile,[]);
    title(['Target Excitation Profile'])
end

% For a more sophisticated Gaussian curve use the following code:
% a = (cos(theta)/sigma_x)^2 + (sin(theta)/sigma_y)^2;
% b = -sin(2*theta)/(sigma_x)^2 + sin(2*theta)/(sigma_y)^2 ;
% c = (sin(theta)/sigma_x)^2 + (cos(theta)/sigma_y)^2;
% [X, Y] = meshgrid(-5:.1:5, -5:.1:5);
% Z = A*exp( - (a*(X-x0).^2 + b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
