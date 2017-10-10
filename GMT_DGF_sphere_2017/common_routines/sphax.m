function sphax(radius,gridflag,axlab,lc,lw)

% plots great circles and cartesian axes on a sphere.
%
% DKS 5/18/93
%
% input arguments:
%       radius: radius of sphere
%       gridflag: 0 --> plot only cartesian axes and great circles (default).
%                 1 --> plot spherical grid.
%       axlab:  0 --> no axis labels
%               1 --> label x,y,z axes
%
% The axes and viewpoint of the plot may be controlled with the Matlab4.0
% commands AXIS and VIEW.
%

drfac=pi/180;

%fill out argument list
if nargin < 5,
    lw = 2;
    if nargin < 4,
        lc = 'k';
        if nargin < 3,
            axlab = 1;
            if nargin < 2,
                gridflag = 0;
            end
        end
    end
end

%plot great circles
ang=2*pi*(0:180)'/180;
%xy
xc=radius*sin(ang); yc=radius*cos(ang); zc=radius*zeros(181,1);
plot3(xc,yc,zc,lc)
axis([-radius radius -radius radius -radius radius])
axis('image')
hold on
%xz
plot3(xc,zc,yc,lc)
%yz
plot3(zc,xc,yc,lc)
if axlab,
        xlabel('x')
        ylabel('y')
        zlabel('z')
end

%plot spherical grid
if gridflag,
        nphi=12;
        for iphi=1:nphi,
                phi=2*pi*(iphi-1)/nphi;
                plot3(xc*cos(phi),xc*sin(phi),yc,':')
        end
        ntheta=6;
        for itheta=1:ntheta,
                theta=pi*(itheta-1)/ntheta;
                plot3(yc*sin(theta),xc*sin(theta),cos(theta)*ones(181,1),':')
        end
end

%plot axes
plot3(radius*[-1 -1/180 1/180 1]',zeros(4,1),zeros(4,1),lc)
plot3(zeros(4,1), radius*[-1 -1/180 1/180 1]', zeros(4,1),lc)
plot3(zeros(4,1),zeros(4,1),radius*[-1 -1/180 1/180 1]',lc)

hold off
