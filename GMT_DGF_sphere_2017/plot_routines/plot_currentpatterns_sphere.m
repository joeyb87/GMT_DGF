function    plot_currentpatterns_sphere(current_id,fieldstrength,plot_opts,path_opts,currentpattern,currentphi,currenttheta,plotlabel,tickint)


scale_axes = 1;
currentpatternmatrixsize = size(currentphi);

figure
if plot_opts.real_part_flag,
    quiver(real(currentpattern(:,:,1)),real(currentpattern(:,:,2)),'LineWidth',plot_opts.currentarrowswidth,'Color',plot_opts.currentarrowscolor); % x-axis = theta, y-axis = phi
    if current_id == 1
        title('Ideal Current Patterns (Real Part)');
    else
        title('Coil Current Patterns (Real Part)');
    end
else
    quiver(imag(currentpattern(:,:,1)),imag(currentpattern(:,:,2)),'LineWidth',plot_opts.currentarrowswidth,'Color',plot_opts.currentarrowscolor); % x-axis = theta, y-axis = phi
    if current_id == 1
        title('Ideal Current Patterns (Imaginary Part)');
    else
        title('Coil Current Patterns (Imaginary Part)');
    end
end
axis tight
delta_phi = (max(currentphi(:))-min(currentphi(:)))/currentpatternmatrixsize(1);
delta_theta = (max(currenttheta(:))-min(currenttheta(:)))/currentpatternmatrixsize(2);
if scale_axes,
    set(gca,'DataAspectRatio',[delta_phi delta_theta 1]);
end
rdfac = 180/pi;

deltax = currentpatternmatrixsize(1)/tickint(2); % 64/8 = 8
xtickset = round(1:deltax:currentpatternmatrixsize(1));
xticklab = rdfac*currenttheta(1,xtickset);

deltay = currentpatternmatrixsize(2)/tickint(1); % 64/4 = 16
ytickset = round(1:deltay:currentpatternmatrixsize(2));
yticklab = rdfac*currentphi(ytickset,1);

set(gca,'xtick',xtickset,'ytick',ytickset,...
    'xticklabel',num2str(xticklab(:),3),...
    'yticklabel',num2str(yticklab(:),3))
xlabel('\theta (deg)')
ylabel('\phi (deg)')
set(gcf,'name',plotlabel);


if plot_opts.plot_current_3D,    
    % currentpattern(:,:,1) is the theta component
    if plot_opts.real_part_flag,
        current_x = cos(currentphi).*cos(currenttheta).*real(currentpattern(:,:,1)) - sin(currentphi).*real(currentpattern(:,:,2));
        current_y = sin(currentphi).*cos(currenttheta).*real(currentpattern(:,:,1)) + cos(currentphi).*real(currentpattern(:,:,2));
        current_z = -sin(currenttheta).*real(currentpattern(:,:,1));
    else
        current_x = cos(currentphi).*cos(currenttheta).*imag(currentpattern(:,:,1)) - sin(currentphi).*imag(currentpattern(:,:,2));
        current_y = sin(currentphi).*cos(currenttheta).*imag(currentpattern(:,:,1)) + cos(currentphi).*imag(currentpattern(:,:,2));
        current_z = -sin(currenttheta).*imag(currentpattern(:,:,1));
    end
    %--- OLD rotation matrix ---%
    % current_x = -sin(currentphi).*real(currentpattern(:,:,2)) + cos(currentphi).*cos(currenttheta).*real(currentpattern(:,:,1));
    % current_y = cos(currentphi).*real(currentpattern(:,:,2)) + sin(currentphi).*cos(currenttheta).*real(currentpattern(:,:,1));
    % current_z = -sin(currenttheta).*real(currentpattern(:,:,1));
    %---------------------------%

    [x,y,z] = sphere(size(currentphi,1) - 1);

    figure;
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
        
        %         textfs = 26;                        % for plots used in ISMRM talk
        %         textlc = 'w';                       % for plots used in ISMRM talk
        %         axislc = 'y';                       % for plots used in ISMRM talk
        %         axislw = 2.5;                       % for plots used in ISMRM talk
        %         set(s,'facelighting','phong');      % for plots used in ISMRM talk
        %         set(gcf,'color','k');               % for plots used in ISMRM talk
        %
        %         hlx = line([0 plotlim],[0 0],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
        %         htx = text(1.07*plotlim,0,0,'x','fontsize',textfs,'color',textlc);
        %         hly = line([0 0],[0 plotlim],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
        %         hty = text(0,1.07*plotlim,0,'y','fontsize',textfs,'color',textlc);
        %         hlz = line([0 0],[0 0],[0 plotlim],'linewidth',axislw,'color',axislc,'erasemode','normal');
        %         htz = text(0,0,1.07*plotlim,'z','fontsize',textfs,'color',textlc);
        
        
        if current_id == 1
            
            hlx = line([0 plotlim],[0 0],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
            htx = text(1.07*plotlim,0,0,'x','fontsize',textfs,'color',textlc);
            hly = line([0 0],[0 plotlim],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
            hty = text(0,1.07*plotlim,0,'y','fontsize',textfs,'color',textlc);
            hlz = line([0 0],[0 0],[0 plotlim],'linewidth',axislw,'color',axislc,'erasemode','normal');
            htz = text(0,0,1.07*plotlim,'z','fontsize',textfs,'color',textlc);
            
            view(plot_opts.sphereview)
        else % RL 9/26/2013 artificially fix the 90 degree shift between FOV seen by coil and ultimate calculations
            
            
            hlx = line([0 plotlim],[0 0],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
            htx = text(1.07*plotlim,0,0,'y','fontsize',textfs,'color',textlc);
            hly = line([0 0],[-plotlim 0],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
            hty = text(0,-1.07*plotlim,0,'x','fontsize',textfs,'color',textlc);
            hlz = line([0 0],[0 0],[0 plotlim],'linewidth',axislw,'color',axislc,'erasemode','normal');
            htz = text(0,0,1.07*plotlim,'z','fontsize',textfs,'color',textlc);
            
            viewopts = plot_opts.sphereview;
            viewopts(1) = viewopts(1) - 90;
            view(viewopts)
        end
        
        %         view(plot_opts.sphereview)
        
        
        %         hlx = line([0 -plotlim],[0 0],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
        %         htx = text(-1.07*plotlim,0,0,'x','fontsize',textfs,'color',textlc);
        %         hly = line([0 0],[0 plotlim],[0 0],'linewidth',axislw,'color',axislc,'erasemode','normal');
        %         hty = text(0,1.07*plotlim,0,'y','fontsize',textfs,'color',textlc);
        %         hlz = line([0 0],[0 0],[0 -plotlim],'linewidth',axislw,'color',axislc,'erasemode','normal');
        %         htz = text(0,0,-1.07*plotlim,'z','fontsize',textfs,'color',textlc);
        %         view(-37.5-180,30)
    end
    
    hold off
    if plot_opts.real_part_flag,
        title('Ideal Current Patterns (Real Part)');
    else
        title('Ideal Current Patterns (Imaginary Part)');
    end
    set(gcf,'name',plotlabel);
end
