function plot_polygon(poly, color)

x = poly(:, 1);
y = poly(:, 2);

d = size(poly, 2);
if d == 2
    plot([x; x(1)], [y; y(1)], color);
    
    for i=1:length(x)
       text(x(i), y(i), num2str(i)); 
    end
else
    z = poly(:, 3);
    plot3([x; x(1)], [y; y(1)], [z; z(1)], color);
    hold on;
    plot3([x; x(1)], [y; y(1)], [z; z(1)], '.');
    
    for i=1:length(x)
       text(x(i), y(i), z(i), num2str(i)); 
    end
end


    
    
axis equal;
