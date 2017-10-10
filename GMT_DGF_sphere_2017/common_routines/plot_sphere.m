function plot_sphere(r, theta, phi,c)


x=r*sin(theta)*cos(phi);
y=r*sin(theta)*sin(phi);
z=r*cos(theta);
plot3(x,y,z,c);
hold on;
end