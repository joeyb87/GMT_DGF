function [] = plot_3planes(vec,xcut,ycut,zcut,windowlims)

figure

% -------------------------------------------------------------------------
% Subplot 1

subplot(1,3,1);
imshow(squeeze(vec(xcut,:,:)),windowlims);
colormap('hot');
colorbar;
axis image

% -------------------------------------------------------------------------
% Subplot 2

subplot(1,3,2);
imshow(squeeze(vec(:,ycut,:)),windowlims);
colormap('hot');
colorbar;
axis image


% -------------------------------------------------------------------------
% Subplot 3

subplot(1,3,3);
imshow(squeeze(vec(:,:,zcut)),windowlims);
colormap('hot');
colorbar;
axis image
