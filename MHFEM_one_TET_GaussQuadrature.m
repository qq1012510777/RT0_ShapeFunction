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
    Area_eachFace(j) = A_j;
end

V_ref = 0.5 * 1 / 3;
FaceArea_ref = [Area_tri([1, 0, 0], [0, 1, 0], [0, 0, 1])
    Area_tri([0, 0, 0], [0, 1, 0], [0, 0, 1])
    Area_tri([0, 0, 0], [1, 0, 0], [0, 0, 1])
    Area_tri([0, 0, 0], [1, 0, 0], [0, 1, 0])];

[xa,ya,za,wt]=TetQuadDat(1);
NumQuaPnts = size(xa, 2);

J = [coord(2, 1) - coord(1, 1), coord(3, 1) - coord(1, 1), coord(4, 1) - coord(1, 1);
    coord(2, 2) - coord(1, 2), coord(3, 2) - coord(1, 2), coord(4, 2) - coord(1, 2);
    coord(2, 3) - coord(1, 3), coord(3, 3) - coord(1, 3), coord(4, 3) - coord(1, 3);];
detJ = abs(det(J));

K_permea = eye(3);

for l = 1:NumQuaPnts
    for m = 1:NumQuaPnts
        xi = xa(l);
        eta = ya(l);
        zeta = za(l);
        
        Psi_ref_l = 1 / detJ .* inv(K_permea) * J * [FaceArea_ref(1) / (3 * V_ref) .* ([xi, eta, zeta] - [0, 0, 0])
            FaceArea_ref(2) / (3 * V_ref) .* ([xi, eta, zeta] - [1, 0, 0])
            FaceArea_ref(3) / (3 * V_ref) .* ([xi, eta, zeta] - [0, 1, 0])
            FaceArea_ref(4) / (3 * V_ref) .* ([xi, eta, zeta] - [0, 0, 1])]';
        % Psi_ref_l = 1 / detJ .* inv(K_permea) * J * [3^0.5 .* ([xi, eta, zeta] - [0, 0, 0])
        %      ([-1 + xi, eta, zeta] - [0, 0, 0])
        %      ([xi, -1 + eta, zeta] - [0, 0., 0])
        %      ([xi, eta, -1 + zeta] - [0, 0, 0])]';

        xi_m = xa(m);
        eta_m = ya(m);
        zeta_m = za(m);
        
        Psi_ref_m = [FaceArea_ref(1) / (3 * V_ref) .* ([xi_m, eta_m, zeta_m] - [0, 0, 0])
            FaceArea_ref(2) / (3 * V_ref) .* ([xi_m, eta_m, zeta_m] - [1, 0, 0])
            FaceArea_ref(3) / (3 * V_ref) .* ([xi_m, eta_m, zeta_m] - [0, 1, 0])
            FaceArea_ref(4) / (3 * V_ref) .* ([xi_m, eta_m, zeta_m] - [0, 0, 1])]';
        % Psi_ref_m = 1 / detJ .* J * [3^0.5 .* ([xi_m, eta_m, zeta_m] - [0, 0, 0])
        %      ([-1 + xi_m, eta_m, zeta_m] - [0, 0, 0])
        %      ([xi_m, -1 + eta_m, zeta_m] - [0, 0., 0])
        %      ([xi_m, eta_m, -1 + zeta_m] - [0, 0, 0])]';


        B = B + wt(l) * wt(m) * Psi_ref_l' * Psi_ref_m;
    end
end

C = zeros(4, 1);
b = zeros(5, 1);
% for l = 1:NumQuaPnts
%     xi = xa(l);
%     eta = ya(l);
%     zeta = za(l);
%     C = C + wt(l) * 1/detJ .* [3 * 3^0.5;3;3;3];
% 
%     b(2) = b(2) -wt(l) * 1/detJ .* FaceArea_ref(2) * 100;
%     b(3) = b(3) -wt(l) * 1/detJ .* FaceArea_ref(3) * 1;
% 
% end

C=Area_eachFace;
b(2) = -Area_eachFace(2) * 100;
b(3) = -Area_eachFace(3) * 1;

K = [B, -C;
    -C', 0];


K(4, :) = 0;
K(:, 4) = 0;
K(4, 4) = 1;
b(4, 1) = -0;

K(1, :) = 0;
K(:, 1) = 0;
K(1, 1) = 1;
b(1, 1) = -0;

x = inv(K ) * b
x(1:4) .* Area_eachFace
    % 1.0000         0         0         0         0
    %      0    0.0010         0         0   -0.0555
    %      0         0    0.0007         0   -0.0786
    %      0         0         0    1.0000         0
    %      0   -0.0555   -0.0786         0         0


   %  1.0e+03 *
   % 
   %       0
   % -4.0099
   %  2.8354
   %       0
   %  0.0258