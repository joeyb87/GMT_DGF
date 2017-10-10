function    movie_currentpatterns_sphere(fieldstrength,plot_opts,g_opts,path_opts,currentpattern,currentphi,currenttheta,plotlabel)

% currentpatternmatrixsize = size(currentphi);

delta_time = 1/plot_opts.movie_frames;
time_vect = (0:delta_time:1)./(fieldstrength*42.576e6);

% time_vect = (0:0.01:1)./(fieldstrength*42.576e6);
%   time_vect = 0:100;
omega_value = 2*pi*fieldstrength*42.576e6;
exp_time_term = exp(1i*omega_value*time_vect);

[x,y,z] = sphere(size(currentphi,1) - 1);

h0 = warndlg('Don''t interact with the computer while the movie frames are recorded','** WARNING **');
if 1
set( h0, 'Visible', 'off' );
fontName = 'FixedWidth';
fontSize = 16;
% get handles to the UIControls ([OK] PushButton) and Text
kids0 = findobj( h0, 'Type', 'UIControl' );
kids1 = findobj( h0, 'Type', 'Text' );
% change the font and fontsize
extent0 = get( kids1, 'Extent' );       % text extent in old font
set( [kids0, kids1], 'FontName', fontName, 'FontSize', fontSize );
extent1 = get( kids1, 'Extent' );       % text extent in new font
% need to resize the msgbox object to accommodate new FontName 
% and FontSize
delta = extent1 - extent0;              % change in extent
pos = get( h0, 'Position' );     % msgbox current position
pos = pos + delta;                      % change size of msgbox
set( h0, 'Position', pos );      % set new position
set( h0, 'Visible', 'on' );
end
uiwait(h0);

figure;
pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) pos(3)*plot_opts.movie_window_scaling pos(4)*plot_opts.movie_window_scaling]);

for iframe = 1:length(time_vect)
    net_current = currentpattern*exp_time_term(iframe);
    
    if plot_opts.real_part_flag,
        current_x = cos(currentphi).*cos(currenttheta).*real(net_current(:,:,1)) - sin(currentphi).*real(net_current(:,:,2));
        current_y = sin(currentphi).*cos(currenttheta).*real(net_current(:,:,1)) + cos(currentphi).*real(net_current(:,:,2));
        current_z = -sin(currenttheta).*real(net_current(:,:,1));
    else
        current_x = cos(currentphi).*cos(currenttheta).*imag(net_current(:,:,1)) - sin(currentphi).*imag(net_current(:,:,2));
        current_y = sin(currentphi).*cos(currenttheta).*imag(net_current(:,:,1)) + cos(currentphi).*imag(net_current(:,:,2));
        current_z = -sin(currenttheta).*imag(net_current(:,:,1));
    end
    
    quiver3(-x',-y',-z',-current_x,-current_y,-current_z,3,'LineWidth',plot_opts.currentarrowswidth,'Color',plot_opts.currentarrowscolor); % BW optimized sphere
    hold on;
    s = surf(x,y,z);
    set(s,'erasemode','normal')
    set(s,'facecolor','interp')
    set(s,'FaceAlpha',plot_opts.spherefacealpha);
    set(s,'FaceColor',plot_opts.spherefacecolor);
    set(s,'EdgeColor',plot_opts.sphereedgecolor);
    axis equal
    axis tight
    axis off
    set(gcf,'Color',plot_opts.plot_bkgcolor);
    
    if plot_opts.show_axes
        plotlim = 1.4;
        axislc = plot_opts.axes_color;
        axislw = 2;
        textfs = 14;
        textlc = 'k';
        
        hlx = line([0 plotlim],[0 0],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
        htx = text(1.07*plotlim,0,0,'x','fontsize',textfs,'color',textlc);
        hly = line([0 0],[0 plotlim],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
        hty = text(0,1.07*plotlim,0,'y','fontsize',textfs,'color',textlc);
        hlz = line([0 0],[0 0],[0 plotlim],'linewidth',axislw,'color',axislc,'erasemode','normal');
        htz = text(0,0,1.07*plotlim,'z','fontsize',textfs,'color',textlc);
        
    end
    
    view(plot_opts.sphereview)
    
    disp(['Frame = ' num2str(iframe) ]);
    F(iframe) = getframe(gcf); % [left bottom width height]
    hold off
end
disp('Start making movie...');
counter = 0;
for frame = 1:1:size(F,2)
    counter = counter + 1;
    [X(:,:,:,counter),Map] = frame2im(F(frame));
end
mov = immovie(X,Map);
moviefilename = [path_opts.moviedir '/3D_' plotlabel '.avi'];
movie2avi(mov,moviefilename,'compression',plot_opts.movie_compression,'fps',plot_opts.movie_fps);
disp('movie saved.');

