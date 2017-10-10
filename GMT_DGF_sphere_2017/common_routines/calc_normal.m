function normal = calc_normal(poynting_radius, avg_theta, avg_phi)
x_s=poynting_radius*cos(avg_phi)*sin(avg_theta);
y_s=poynting_radius*sin(avg_phi)*sin(avg_theta);
z_s=poynting_radius*cos(avg_theta);
normal_denominator=sqrt(x_s^2+y_s^2+z_s^2);
normal=1/normal_denominator*[x_s y_s z_s];
%         normal=[normal_1;normal_1];
%         quiver3(x_s,y_s,z_s,normal(1), normal(2), normal(3))
end % end calc_normal
