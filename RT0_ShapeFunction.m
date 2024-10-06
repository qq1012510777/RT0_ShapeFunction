clc
clear all
close all

% coordinates of the source element
coord = [    3.2217    9.8066
    4.6208    8.1898
    4.9911   10.3278];

figure(1)
subplot(1, 2, 1)
title("Source triangle")
patch('Vertices', coord, 'Faces', 1:3, 'FaceVertexCData', 0, 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
xlabel('x')
ylabel('y')
pbaspect([1, 1, 1])
for i = 1:3
    text(coord(i, 1), coord(i, 2), ['(', num2str(i), ')'], 'color', 'r'); hold on
end

% the basis functions of RT0 for reference elements is
% \hat{\psi} = a + b \hat{x}
% for the reference triangular element, the coordinates are
coord_ref = [0, 0;
    1, 0;
    0, 1];
figure(1)
subplot(1, 2, 2)
title("Reference triangle")
patch('Vertices', coord_ref, 'Faces', 1:3, 'FaceVertexCData', 0, 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
xlabel('x')
ylabel('y')
pbaspect([1, 1, 1])
for i = 1:3
    text(coord_ref(i, 1), coord_ref(i, 2), ['(', num2str(i), ')'], 'color', 'r'); hold on
end

% we have to calculate the values of a and b for each face (or say, edge)
% of the reference element. Note that face 1 is the edge opposite node 1
% we use the formula $\frac{1}{|f_i| \int_{f_i} \hat{psi_i}(x_j) \cdot \n_j \mathrm d s = \delta_{ij}$
% where $f_i$ is the $i$th face, and $|f_i|$ is the length of the face,
% $x_j$ here is the center of the face, $n_j$ is the outer normal vector on
% the face, and $\delta _{ij}$ is the Kronecker delta.

% now let us calculate a and b for shape functions on reference elements
% first, the centers of face 1, 2, and 3
CenterFaces_ref = [(coord_ref(2, :) + coord_ref(3, :)) .* 0.5
    (coord_ref(1, :) + coord_ref(3, :)) .* 0.5
    (coord_ref(1, :) + coord_ref(2, :)) .* 0.5];
for i = 1:3
    subplot(1, 2, 2)
    text(CenterFaces_ref(i, 1), CenterFaces_ref(i, 2), [num2str(i)], 'color', 'b'); hold on
end

% the normal outer vector of each face
% and the face length
OuterNormalFaces_ref = [coord_ref(3, :) - coord_ref(2, :);
    coord_ref(1, :) - coord_ref(3, :);
    coord_ref(2, :) - coord_ref(1, :)];
FaceLength_ref = [vecnorm(OuterNormalFaces_ref')]';
OuterNormalFaces_ref = [OuterNormalFaces_ref(:, 2), -OuterNormalFaces_ref(:, 1)];
OuterNormalFaces_ref = OuterNormalFaces_ref ./ [vecnorm(OuterNormalFaces_ref')]';
subplot(1, 2, 2)
quiver(CenterFaces_ref(:, 1), CenterFaces_ref(:, 2), OuterNormalFaces_ref(:, 1), OuterNormalFaces_ref(:, 2), 'color', 'g'); hold on

% calculate a and b
% inv([1/2^0.5, 1/2^0.5, 1/2^0.5; -1, 0, 0; 0, -1, 0]) * [1, 0, 0; 0, 1, 0; 0, 0, 1]
syms a1 a2 b n1 n2 real
coe_matrix = zeros(3, 3);
for i = 1:3
    dot_product = dot([a1 + b * CenterFaces_ref(i, 1), a2 + b * CenterFaces_ref(i, 2)], [n1, n2]);
    dot_product = collect(dot_product, [a1, a2, b]);
    Ax = coeffs(dot_product, [b, a2, a1]);
    coe_matrix(i, :) = subs(Ax, [n1, n2], OuterNormalFaces_ref(i, :));
end
coe_matrix  = inv(coe_matrix) * eye(3);

for i = 1:3
    disp(['psi_hat_', num2str(1), ' = [', num2str(coe_matrix(1, i)), ' + ', num2str(coe_matrix(3, i)), ' * x_hat', '; ' num2str(coe_matrix(2, i)), ' + ', num2str(coe_matrix(3, i)), ' * y_hat', ']'])
end

% Piola transformation
% calculate $\psi_i$ of source elements
B_K = [coord(2, 1) - coord(1, 1), coord(3, 1) - coord(1, 1);
    coord(2, 2) - coord(1, 2), coord(3, 2) - coord(1, 2)];
% validate B_K
disp(['validate B_K, see if the result is [0;0]'])
for i = 1:3
    % see if the result is [0;0]
    disp(B_K * coord_ref(i, :)' + coord(1, :)' - coord(i, :)');
end

detB_k = det(B_K);
Area_T = Area_tri(coord(1, :), coord(2, :), coord(3, :));

CloudPnts = GenCloudPointsInTriangle(coord(1, :), coord(2, :), coord(3, :), 40);
% subplot(1, 2, 1)
% scatter(CloudPnts(:, 1), CloudPnts(:, 2), '+');
% hold on

q_n = [-2 / norm(coord(2, :) - coord(3, :)), 1 / norm(coord(2, :) - coord(3, :)), 1/ norm(coord(3, :) - coord(1, :))];
V_t = zeros(size(CloudPnts, 1), 2);
for i = 1:size(CloudPnts, 1)
    coord_hat = inv(B_K) * [CloudPnts(i, :) - coord(1, :)]';
    xi_s = coord_hat(1);
    eta_s = coord_hat(2);
    
    M_j = zeros(2, 1);

    for j = 1:3
        psi_hat_xi = coe_matrix(1, j) + coe_matrix(3, j) * xi_s;
        psi_hat_eta = coe_matrix(2, j) + coe_matrix(3, j) * eta_s;

        M_j = M_j + 1 / detB_k .* B_K * [psi_hat_xi, psi_hat_eta]' .* q_n(j);
    end
    V_t(i, :) = M_j';

end

subplot(1, 2, 1)
quiver(CloudPnts(:, 1) , CloudPnts(:, 2), V_t(:, 1), V_t(:, 2)); hold on

B = zeros(3, 3);
C = zeros(3, 1);

xi = [1/2., 1/2., 0];
eta = [1/2., 0, 1/2.];
weights = [1/6., 1/6., 1/6.];

for i = 1:3
    for j = 1:3
        for l = 1:3
            for m = 1:3
                psi_hat_i_xi = coe_matrix(1, i) + coe_matrix(3, i) * xi(l);
                psi_hat_i_eta = coe_matrix(2, i) + coe_matrix(3, i) * eta(l);
                
                psi_hat_j_xi = coe_matrix(1, j) + coe_matrix(3, j) * xi(m);
                psi_hat_j_eta = coe_matrix(2, j) + coe_matrix(3, j) * eta(m);

                w1 = weights(l);
                w2 = weights(m);

                B(i, j) = B(i, j) + (1 / detB_k)^2 * w1 * w2 * dot(B_K * [psi_hat_i_xi, psi_hat_i_eta]', B_K * [psi_hat_j_xi, psi_hat_j_eta]');
            end
            %--------------------------------
        end
    end
end

div_phi_hat = [2^0.5 + 2^0.5; 1 + 1; 1 + 1];

 
C = (-1 / detB_k) * div_phi_hat .* Area_tri(coord_ref(1, :), coord_ref(2, :), coord_ref(3, :));


B
C
%1.04994164052243	-0.260329062344791	-0.0499426798279379
%-0.260329062344791	0.965009678007263	-0.281542085802101
%-0.0499426798279379	-0.281542085802101	1.04847072628005

K = [B, C;C', 0];
b = [norm(coord(2, :) - coord(3, :)) * 100;
    norm(coord(1, :) - coord(3, :)) * 10;
    norm(coord(1, :) - coord(2, :)) * 10;
    0];

inv(K) * b