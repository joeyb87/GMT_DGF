function [] = plot_fields(field3d,x_fov,y_fov,z_fov,mask_reg1,mask_reg2,mask_reg3,mask_reg4,s_opts)

[Nx,Ny,Nz] = size(x_fov);
scat = struct;
scat.mask = mask_reg4;
scat.ind = find(scat.mask);
scat.Ns = numel(scat.ind);
scat.dom.res = 1;
scat.dom.dims = [Nx,Ny,Nz];
off = ceil(scat.dom.dims/2);
x = scat.dom.res*((1:Nx).'-off(1));
y = scat.dom.res*((1:Ny).'-off(2));
z = scat.dom.res*((1:Nz).'-off(3));
[rX,rY,rZ] = ndgrid(x,y,z);
scat.dom.r = reshape([rX(:),rY(:),rZ(:)],[Nx,Ny,Nz,3]);
% maps = {map1,map2;map3,[]};
% titles = {'title1','title2','title3',''};

netfield3d = squeeze(sum(field3d,1));
mag_netfield3d = sqrt( abs(netfield3d(:,1)).^2 + abs(netfield3d(:,2)).^2 + abs(netfield3d(:,3)).^2 );
mag_netfield3d = reshape(mag_netfield3d,[Nx,Ny,Nz]);
maps = {mag_netfield3d(scat.ind),[]};
titles = {'B basis',''};
% 
% lims = {[0,1],[0,1];[0,1],[0,1]};
% f = slice_view(scat,maps,titles,lims);
f = slice_view(scat,maps,titles);