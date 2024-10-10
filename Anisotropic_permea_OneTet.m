clc
clear all
close all

coord = [         0    0.3333         0
    0.3333    0.3333    0.3333
         0    0.3333    0.3333
         0         0         0];
element = [1, 2, 3, 4];

faceNode = [2, 3, 4; 3, 4, 1; 4, 1, 2; 1, 2, 3];

V = tetrahedronVolume(coord(1, :), coord(2, :), coord(3, :), coord(4, :));

B = zeros(4, 4);
C = zeros(4, 1);

K = eye(3, 3);
delta_lm = eye(4);

for i = 1:4
    NodeFace_i = faceNode(i, :);
    A_i = Area_tri(coord(NodeFace_i(1), :), coord(NodeFace_i(2), :), coord(NodeFace_i(3), :));

    for j = 1:4
        NodeFace_j = faceNode(j, :);
        A_j = Area_tri(coord(NodeFace_j(1), :), coord(NodeFace_j(2), :), coord(NodeFace_j(3), :));

        for k = 1:3       
            for l = 1:3  

                for a = 1:4
                    for b = 1:4
                        B(i, j) = B(i, j) + A_i * A_j / (9 * V^2) * ...
                            K(k, l) * (coord(a, k) - coord(i, k)) * (coord(b, l) - coord(j, l)) * V / 20 .* (1 + 1 * delta_lm(a, b));
                    end
                    
                end
            end
        end
    end
end
