close all;

n = 50;

hi = 1/n;
hj = 1/n;

I=1:n;
J=1:n;

Ni = length(I);
Nj = length(J);

[X, Y] = meshgrid(I*hi, J*hj);

idx = @(i, j) (i-1) * Nj + j;

ndofs = Ni * Nj;

A = sparse(ndofs, ndofs);

h2  = hi*hj;
hi2 = hi*hi;
hj2 = hj*hj;

is_b_dof = zeros(ndofs, 1);

for i=I
    for j=J
        k = idx(i, j);
        A(k, k) = 4.0/h2;
        
        is_boundary = 0;
        
        if i+1 <= Ni
           l = idx(i + 1, j);
           A(k, l) = -1/hi2;
        else
            is_boundary = 1;
        end
        
        if i-1 > 0
            l = idx(i - 1,j);
            A(k, l) = -1/hi2;
        else
            is_boundary = 1;
        end
        
        if j+1 <= Nj
           l = idx(i, j + 1);
           A(k, l) = -1/hj2;
        else
            is_boundary = 1;
        end
        
        if j-1 > 0
            l = idx(i, j - 1);
            A(k, l) = -1/hj2;
        else
            is_boundary = 1;
        end
        
        is_b_dof(k) = is_boundary;
    end 
end

A = A*A*A;

for k=1:ndofs
  if is_b_dof(k) == 1
     A(k, :) = 0; 
     A(k, k) = 1;
  end
end

% rhs = ones(ndofs, 1);
% rhs(bidx) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exponential function

C = ones(size(X))*0.5;
DX = X - C;
DY = Y - C;

D = sqrt(DX.*DX + DY.*DY);

alpha = 10;
beta = 50;
E = alpha * exp(-beta * D);
figure;
mesh(X, Y, E);

z = E(:);
rhs = z;
rhs(bidx) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non-linear function 
% Z = sin(30 * X.*X + Y.*Y);
% z = Z(:);
% laplz = A*z;
% rhs = laplz;
% rhs(bidx) = z(bidx);
% 
% figure;
% mesh(X, Y, Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve

u = A\rhs;
U = reshape(u, size(Z));

figure;
mesh(X, Y, U);

% isB = reshape(is_b_dof, size(Z));
% 
% figure;
% mesh(X, Y, isB);
