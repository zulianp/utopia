function v = tet10(i, p)
    switch(i)
        case 0 
            v = f0(p);
        case 1 
            v = f1(p);
        case 2 
            v = f2(p);
        case 3 
            v = f3(p);
        case 4 
            v = f4(p);
        case 5 
            v = f5(p);
        case 6 
            v = f6(p);
        case 7 
            v = f7(p);
        case 8 
            v = f8(p);
        case 9 
            v = f9(p);
    end
return

function [v] = f0(p)
L0 = 1 - p(1, :) - p(2, :) - p(3, :);
v = (2.0 .* L0 - 1.0) .* L0;
return


function [v] = f1(p)
L1 = p(1, :);
v = (2.0 .* L1 - 1.0) .* L1;
return


function [v] = f2(p)
L2 = p(2, :);
v = (2.0 .* L2 - 1.0) .* L2;
return


function [v] = f3(p)
L3 = p(3, :);
v = (2.0 .* L3 - 1.0) .* L3;
return

function [v] = f4(p)
L0 = 1 - p(1, :) - p(2, :) - p(3, :);
L1 = p(1, :);
v = 4.0 .* L0 .* L1;
return


function [v] = f5(p)
L1 = p(1, :);
L2 = p(2, :);
v = 4.0 .* L1 .* L2;
return


function [v] = f6(p)
L0 = 1 - p(1, :) - p(2, :) - p(3, :);
L2 = p(2, :);
v = 4.0 .* L0 .* L2;
return


function [v] = f7(p)
L0 = 1 - p(1, :) - p(2, :) - p(3, :);
L3 = p(3, :);
v = 4.0 .* L0 .* L3;
return


function [v] = f8(p)
L1 = p(1, :);
L3 = p(3, :);
v = 4.0 .* L1 .* L3;
return


function [v] = f9(p)
L2 = p(2, :);
L3 = p(3, :);
v = 4.0 .* L2 .* L3;
return
