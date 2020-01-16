function [dt] = shear_zone_grid(path)
d = load(path);
d = d';
d = d(:);

ntot = length(d);
n = ntot/3;

x = d(1:n);
y = d((n+1):(2*n));
z = d((2*n+1):(3*n));

thick_layer = (max(z) - min(z)) * 0.1;


dt = delaunayTriangulation(x, y);
tri = dt(:,:);

close all; 
trimesh(tri, x, y, z); 

TR = triangulation(tri, [x'; y'; z']');
stlwrite(TR,'shear-zone.stl','text');

boundaryEdges = dt.freeBoundary();
xb = x(boundaryEdges);
yb = y(boundaryEdges);
zb = z(boundaryEdges);
clf;

hold on;
plot3(xb, yb, zb, '-r');

axis equal; axis tight;

%offset triangles
bottom_tri = tri + n;

%swap orientation
bottom_tri = [ 
    bottom_tri(:, 1), ...
    bottom_tri(:, 3), ...
    bottom_tri(:, 2)  ...
]; 

nboundary  = length(boundaryEdges);
side_tri = zeros(2 * nboundary, 3);

for i=1:nboundary
    v1 = boundaryEdges(i, 1);
    v2 = boundaryEdges(i, 2);
    
    idx = (i - 1) * 2 + 1;
    
    side_tri(idx, :)     = [ v1, v2 + n, v2 ];
    side_tri(idx + 1, :) = [ v1, v1 + n, v2 + n];
end

v_bottom_tri  = [tri; bottom_tri; side_tri];
v_bottom_pts  = [x', x'; y', y'; z', (ones(size(z)) * min(z) - thick_layer)']';

TR = triangulation(v_bottom_tri, v_bottom_pts);
stlwrite(TR,'shear-zone-bottom.stl','text');

v_top_pts  = [x', x'; y', y'; (ones(size(z)) * max(z) + thick_layer)', z']';

TR = triangulation(v_bottom_tri, v_top_pts);
stlwrite(TR,'shear-zone-top.stl','text');

end
