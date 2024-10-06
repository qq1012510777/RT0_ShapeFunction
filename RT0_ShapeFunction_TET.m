clc
clear all
close all

coord = [         0    0.3333         0
    0.3333    0.3333    0.3333
         0    0.3333    0.3333
         0         0         0];
element = [1, 2, 3, 4];
faceNode = [2, 3, 4; 3, 4, 1; 4, 1, 2; 1, 2, 3];
figure(1)
view(3)
patch('Vertices', coord, 'Faces', faceNode, ...
         'FaceVertexCData', zeros(size(element, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); 
text(coord(:, 1), coord(:, 2), coord(:, 3), num2str(element'))
hold on
xlabel('x'); ylabel('y'); zlabel('z');
pbaspect([1, 1, 1])
CenterEachFaces = zeros(4, 3);
for i = 1:4
    CenterEachFaces(i, :) = (coord(faceNode(i, 1), :) + coord(faceNode(i, 2), :) + coord(faceNode(i, 3), :)) .* 1.0/3.0;
end
text(CenterEachFaces(:, 1), CenterEachFaces(:, 2), CenterEachFaces(:, 3), num2str([1, 2, 3, 4]'), 'Color', 'g'); hold on

psi_coe = zeros(4, 1);

V_t = tetrahedronVolume(coord(1, :), coord(2, :), coord(3, :), coord(4, :));
Area_each_face = zeros(4, 1);
NormalVector_each_face = zeros(4, 3);
CentersTET = 1/4.0 * (coord(1, :) + coord(2, :) + coord(3, :) + coord(4, :));
for i = 1:4
    Area_t = Area_tri(coord(faceNode(i, ...
         1), :), coord(faceNode(i, ...
         2), :), coord(faceNode(i, ...
         3), :));
    psi_coe(i) = Area_t / (3 * V_t);
    Area_each_face(i) = Area_t;

    AB = coord(faceNode(i, ...
         2), :) - coord(faceNode(i, ...
         1), :);
    AC = coord(faceNode(i, ...
         3), :) - coord(faceNode(i, ...
         2), :);
    NormalVector_each_face(i, :) = cross(AB, AC);
    NormalVector_each_face(i, :) = NormalVector_each_face(i, :) ./ norm(NormalVector_each_face(i, :));

    KO =  CentersTET - CenterEachFaces(i, :);
    KO = KO / norm(KO);
    if (dot(KO, NormalVector_each_face(i, :) ) > 0)
        NormalVector_each_face(i, :) = -NormalVector_each_face(i, :);
    end
end

q_n = [1 / Area_each_face(1), 1 / Area_each_face(2), -1 / Area_each_face(3), -1 / Area_each_face(4)];

points_cloud = generatePointsInTetrahedron(coord(1, :), coord(2, :), coord(3, :), coord(4, :), 40);

figure(1)
scatter3(points_cloud(:, 1), points_cloud(:, 2), points_cloud(:, 3), '+');

V = points_cloud .* 0;
Vcenter = CenterEachFaces .* 0;
for i = 1:4
    V = V + psi_coe(i) .* (points_cloud - coord(i, :)) .* q_n(i);
    Vcenter = Vcenter + psi_coe(i) .* (CenterEachFaces - coord(i, :)) .* q_n(i);
end

figure(1)
quiver3(points_cloud(:, 1), points_cloud(:, 2), points_cloud(:, 3), ...
    V(:, 1), V(:, 2), V(:, 3), 'color', 'b')
quiver3(CenterEachFaces(:, 1), CenterEachFaces(:, 2), CenterEachFaces(:, 3), ...
    Vcenter(:, 1), Vcenter(:, 2), Vcenter(:, 3), 'color', 'r')
quiver3(CenterEachFaces(:, 1), CenterEachFaces(:, 2), CenterEachFaces(:, 3), ...
    NormalVector_each_face(:, 1), NormalVector_each_face(:, 2), NormalVector_each_face(:, 3), 'color', 'c')

for i = 1:4
    V_norm_i = V(i, :);
    V_norm_i = V_norm_i ./ norm(V_norm_i);

    if (dot(V_norm_i, NormalVector_each_face(i, :)) < 0)
        qn = -dot(NormalVector_each_face(i, :), V(i, :)) .* Area_each_face(i)
    else
        qn = -dot(NormalVector_each_face(i, :), V(i, :)) .* Area_each_face(i)
    end

end