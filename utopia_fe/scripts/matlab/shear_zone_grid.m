function [dt] = shear_zone_grid(path)
d = load(path);
d = d';
d = d(:);

ntot = length(d);
n = ntot/3;

x = d(1:n);
y = d((n+1):(2*n));
z = d((2*n+1):(3*n));

dt = delaunayTriangulation(x, y);
tri = dt(:,:);

close all; 
trimesh(tri, x, y, z); 


boundaryEdges = dt.freeBoundary();
xb = x(boundaryEdges);
yb = y(boundaryEdges);
zb = z(boundaryEdges);
clf;

hold on;
plot3(xb, yb, zb, '-r');

axis equal; axis tight;


% fault layer
% TR = triangulation(tri, [x'; y'; z']');
% stlwrite(TR,'shear-zone.stl','text')
% 
% % bottom layer
% bottom = ones(size(z)) * min(z);
% TR = triangulation(tri, [x'; y'; bottom']');
% stlwrite(TR,'flat-zone.stl','text')

end