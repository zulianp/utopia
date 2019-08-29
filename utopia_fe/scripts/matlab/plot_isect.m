function plot_isect(ray, t1, t2)

x = ray(:, 1);
y = ray(:, 2);

xx = [x(1); x(1)];
yy = [y(1); y(1)];
dx = [t1*x(2); t2*x(2)];
dy = [t1*y(2); t2*y(2)];

d =size(ray, 2);
if d == 2
    plot(xx + dx, yy + dy, 'r*');
    quiver(xx, yy, dx, dy);
else
    z = ray(:, 3);
    zz = [z(1); z(1)];
    dz = [t1*z(2); t2*z(2)];

    plot3(xx + dx, yy + dy, zz + dz, 'r*');
    quiver3(xx, yy, zz, dx, dy, dx);
end



axis equal;
len = sqrt(x(2)^2 + y(2)^2);

