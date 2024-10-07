clear all
close all
clc

syms x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4 x y z real

%%%%%% regular triangle
V_t = [x1, x2, x3;
     y1, y2, y3;
     1,  1,  1,];
P_t = [x; y; 1];
lambda_t = V_t \ P_t;
% f = x - x.^2 - x.*y;
% f = 1 - x - y;
f = lambda_t(1) * lambda_t(2);
f = subs(f, [x1, x2, x3, y1, y2, y3], [0, 1, 0, 0, 0, 1]);
result = int(int(f, y, 0, 1 - x), x, 0, 1)

%%%%%% regular tetrahedron
% Define the vertices of the tetrahedron
V = [x1, x2, x3, x4;
     y1, y2, y3, y4;
     z1, z2, z3, z4;
     1,  1,  1,  1];

% Define the point (x, y, z)
P = [x; y; z; 1];

% Solve for barycentric coordinates
lambda = V \ P;
f = lambda_t(1) * lambda_t(3);
f = subs(f, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
integral_result = int(int(int(f, z, 0, 1 - x - y), y, 0, 1 - x), x, 0, 1)