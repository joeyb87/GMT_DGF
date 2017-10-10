function    plot_geometry_sphere(sphere_radius,coil_radii,coil_offsets,coil_rotations,x_fov,y_fov,z_fov,lw,ls,lc,axisflag,labelflag,axislw,axislc,textlc,textfs)

if nargin < 16,
    textfs = 14;
    if nargin < 15,
        textlc = 'k';
        if nargin < 14,
            axislc = 'k';
            if nargin < 13,
                axislw = 2;
                if nargin < 12,
                    labelflag = 1;
                    if nargin < 11,
                        axisflag = 0;
                        if nargin < 10,
                            lc = 'k';
                            if nargin < 9,
                                ls = '-';
                                if nargin < 8,
                                    lw = 2;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% figure; % For g-factor maps subplots
[xs,ys,zs] = sphere(30);
s = surfl(sphere_radius*xs,sphere_radius*ys,sphere_radius*zs);
% s = surfl(0.067*xs,0.067*ys,0.067*zs); % for plots used in ISMRM talk (coils away from sample)

% set(s,'erasemode','normal')
set(s,'facecolor','interp')

% set(s,'facecolor',[0.3 0.3 0.3]) % For g-factor maps subplots

set(s,'facecolor',[0.89 0.89 0.89])    % for BW optimized plots
% set(s,'facecolor',[0 0 0])    % for BW optimized ultimate case plot
set(s,'edgecolor',[0.6 0.6 0.6])    % for BW optimized plots


% set(s,'facecolor',[0.59 0.95 0.97]) % for plots used in ISMRM talk (blue sphere)
% set(s,'edgecolor',[0.28 0.65 1])    % for plots used in ISMRM talk (blue sphere)

% set(s,'facecolor',[1 0.46 0])         % for ultimate case
% set(s,'edgecolor',[0.98 0.62 0.32])   % for ultimate case
% textlc = 'w';                       % for plots used in ISMRM talk
% textfs = 26;                        % for plots used in ISMRM talk
% axislc = 'y';                       % for plots used in ISMRM talk
% axislw = 2.5;                        % for plots used in ISMRM talk

% set(s,'edgecolor','none')
% set(s,'facecolor','none')

% set(s,'edgecolor',[0.83 0.83 0.83])    % for plots used in ISMRM talk
% set(s,'facecolor',[0.41 0.41 0.41])    % for plots used in ISMRM talk
% set(gcf,'color','black');              % for plots used in ISMRM talk
set(s,'facelighting','phong');




plotlim = 1.4*max([sphere_radius max(sqrt(coil_radii.^2 + coil_offsets.^2))]);
axis equal
axis([-plotlim plotlim -plotlim plotlim -plotlim plotlim])
xlabel('x')
ylabel('y')
zlabel('z')
if ~axisflag,
    axis off
end
% title(['\theta=' num2str(coil_rotation(1)) ',\phi=' num2str(coil_rotation(2))])
plotplan_sphere(x_fov,y_fov,z_fov)
ncoils = length(coil_radii);

lc = 'k'; % for BW optimized plots

lc = [1 0.46 0];


for icoil = 1:ncoils,
    [coilvec,nc,ncfull,nseg,combmat] = circular_coil(coil_radii(icoil),coil_offsets(icoil),coil_rotations(icoil,:));
    plotcarr_sphere(coilvec,nc,nseg,lw,ls,lc)
    if labelflag,
        r_coil = 1.1*sphere_radius;
        theta_coil = pi*coil_rotations(icoil,1)/180;
        phi_coil = pi*coil_rotations(icoil,2)/180;
        text(r_coil*sin(theta_coil)*cos(phi_coil),r_coil*sin(theta_coil)*sin(phi_coil),r_coil*cos(theta_coil),num2str(icoil),...
            'fontsize',textfs,'color','y');
    end
end
% hlx = line([0 plotlim],[0 0],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
% htx = text(1.1*plotlim,0,0,'x','fontsize',textfs,'color',textlc);
% hly = line([0 0],[0 plotlim],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
% hty = text(0,1.1*plotlim,0,'y','fontsize',textfs,'color',textlc);
% hlz = line([0 0],[0 0],[0 plotlim],'linewidth',axislw,'color',axislc,'erasemode','normal');
% htz = text(0,0,1.1*plotlim,'z','fontsize',textfs,'color',textlc);

view(-37.5-180,30)
colormap bone
% alpha(.3)
