function plot_segment(poly, color)
x = poly(:, 1);
y = poly(:, 2);

plot(x, y,color);

hold on;

ux = x(2:end) - x(1:end-1);
uy = y(2:end) - y(1:end-1);

bx =  0.5 * (x(2:end) + x(1:end-1));
by =  0.5 * (y(2:end) + y(1:end-1));

lens = sqrt(ux.^2 + uy.^2);
ux = ux ./ lens;
uy = uy ./ lens;

nx = -uy;
ny = ux;

quiver(bx, by, nx, ny);
axis equal;
