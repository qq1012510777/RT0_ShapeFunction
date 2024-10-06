clear all
close all
clc

coord = [         0    0.3333         0
    0.3333    0.3333    0.3333
         0    0.3333    0.3333
         0         0         0];
element = [1, 2, 3, 4];

faceNode = [2, 3, 4; 3, 4, 1; 4, 1, 2; 1, 2, 3];

V = tetrahedronVolume(coord(1, :), coord(2, :), coord(3, :), coord(4, :));

B = zeros(4, 4);
C = zeros(4, 1);

delta_lm = eye(4);

Area_eachFace = zeros(4, 1);

for j = 1:4
    NodeFace_j = faceNode(j, :);
    A_j = Area_tri(coord(NodeFace_j(1), :), coord(NodeFace_j(2), :), coord(NodeFace_j(3), :));
    for k = 1: 4
        NodeFace_k = faceNode(k, :);
        A_k = Area_tri(coord(NodeFace_k(1), :), coord(NodeFace_k(2), :), coord(NodeFace_k(3), :));

        for l = 1:4
            for m = 1:4
                lambda_lm = V / (60 .* (1 + delta_lm(l, m)));
                B(j, k) = B(j, k) + A_j * A_k / (9 * V .^ 2) * dot((coord(l, :) - coord(j, :)), (coord(m, :) - coord(k, :))) * lambda_lm;
            end
        end
        
    end
    C(j) = A_j;
    Area_eachFace(j) = A_j;
end

% in 2
% out 3

b = zeros(5, 1);

b(2) = -Area_eachFace(2) * 20;
b(3) = -Area_eachFace(3) * 1;

K = [B, -C;
    -C', 0];

K(4, :) = 0;
K(:, 4) = 0;
K(4, 4) = 1;
b(4, 1) = -100;

K(1, :) = 0;
K(:, 1) = 0;
K(1, 1) = 1;
b(1, 1) = -100;

inv(K) * b
