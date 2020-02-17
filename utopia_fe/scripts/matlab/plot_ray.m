function plot_ray(ray, color)

x = ray(:, 1);
y = ray(:, 2);

d = size(ray, 2);

if d == 2
    quiver(x(1), y(1), x(2), y(2), 0.1);
else
    z = ray(:, 3);
    quiver3(x(1), y(1), z(1), x(2), y(2), z(2));
end
axis equal;

len = sqrt(x(2)^2 + y(2)^2);



