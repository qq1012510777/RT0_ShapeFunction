clc
clear all
close all

box_tets_3D;

Points = msh.POS;
Element = msh.TETS(:, 1:4);

% Points = [0 0 0;
%     1 0 0;
%     1 1 0;
%     0 1 0;
%     0 0 1;
%     1 0 1;
%     1 1 1;
%     0 1 1;
%     0.5 0.5 0;
%     0.5 0.5 1;
%     0.5 0 0.5;
%     1 0.5 0.5;
%     0.5 1 0.5;
%     0 0.5 0.5;];
% Element = [  9 13 14 12 
%  11 14 10 12 
%  11 9 14 12 
%  14 13 10 12 
%  14 5 1 11 
%  10 5 11 6 
%  8 10 5 14 
%  8 4 13 14 
%  4 1 9 14 
%  13 8 10 7 
%  4 13 9 3 
%  7 6 10 12 
%  12 6 11 2 
%  1 9 11 2 
%  3 7 13 12 
%  12 9 3 2 
%  14 1 9 11 
%  10 5 14 11 
%  10 8 13 14 
%  9 13 4 14 
%  6 11 10 12 
%  13 9 3 12 
%  10 13 7 12 
%  12 11 9 2 ];

figure(1)
subplot(1, 2, 1)
title('Show one tetrahedron')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetNo = 1;
NodeIDforOneTet = [Element(tetNo, [2, 3, 4]);
    Element(tetNo, [3, 4, 1]);
    Element(tetNo, [4, 1, 2]);
    Element(tetNo, [1, 2, 3])];
patch('Vertices', Points, 'Faces', NodeIDforOneTet, ...
         'FaceVertexCData', zeros(size(NodeIDforOneTet, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); 
pbaspect([1, 1, 1])
text(Points(Element(tetNo, :), 1), Points(Element(tetNo, :), 2), Points(Element(tetNo, :), 3), num2str([1;2;3;4]))
CenterEachFace = zeros(4, 3);
for i = 1:4
    CenterEachFace(i, :) = (Points(NodeIDforOneTet(i, 1), :) + Points(NodeIDforOneTet(i, 2), :) + Points(NodeIDforOneTet(i, 3), :)) .* 1.0 / 3.0;
end
text(CenterEachFace(:, 1), CenterEachFace(:, 2), CenterEachFace(:, 3), num2str([1;2;3;4]), 'Color', 'r')

figure(1)
subplot(1, 2, 2)
title('All tetras')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetNo = 1;
NodeIDforOneTet = [Element(tetNo, [2, 3, 4]);
    Element(tetNo, [3, 4, 1]);
    Element(tetNo, [4, 1, 2]);
    Element(tetNo, [1, 2, 3])];
patch('Vertices', Points, 'Faces', [Element(:, [2, 3, 4]); Element(:, [3, 4, 1]); ...
            Element(:, [4, 1, 2]); ...
            Element(:, [1, 2, 3])], ...
         'FaceVertexCData', zeros(size(Element, 1) * 4, 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); 
pbaspect([1, 1, 1])

%%% give global face numbers to each face
FaceLocalOrder = [2, 3, 4
   3, 4, 1
   4, 1, 2
   1, 2, 3];

NumEles = size(Element, 1);
NumPnts = size(Points, 1);

FaceGlobalSort = zeros(NumEles * 4, 3);
NumGlobalFaces = 1;
FaceGlobalID = zeros(NumEles, 4);
FaceAttribute = zeros(NumEles, 4); % 1: x = 0, 2: x = 1, 
                                   % 3: y = 0, 4: y = 1
                                   % 5: z = 0, 6: z = 1
                                   % 0: interior face
x0Face = [];
x1Face = [];
y0Face = [];
y1Face = [];
z0Face = [];
z1Face = [];
for i = 1:NumEles
    for j = 1:4
        P1_ID_j = Element(i, FaceLocalOrder(j, 1));
        P2_ID_j = Element(i, FaceLocalOrder(j, 2));
        P3_ID_j = Element(i, FaceLocalOrder(j, 3));
        
        P = sort([P1_ID_j, P2_ID_j, P3_ID_j]);
        
        rowIndex = find(all(FaceGlobalSort(1:NumGlobalFaces, :) == P, 2));
        if ~isempty(rowIndex)
            if (size(rowIndex, 1) * size(rowIndex, 2) ~= 1)
                fprintf('Wrong\n');
                exit
            end
            FaceGlobalID(i, j) = rowIndex;
        else
            FaceGlobalSort(NumGlobalFaces, :) = P;
            FaceGlobalID(i, j) = NumGlobalFaces;
            NumGlobalFaces = NumGlobalFaces + 1;
        end
        
        if all(Points(P, 1) == 0)
            FaceAttribute(i, j) = 1;
            x0Face = [x0Face ; P];
        elseif all(Points(P, 1) == 1)
            FaceAttribute(i, j) = 2;
            x1Face = [x1Face ; P];
        elseif all(Points(P, 2) == 0)
            FaceAttribute(i, j) = 3;
            y0Face = [y0Face ; P];
        elseif all(Points(P, 2) == 1)
            FaceAttribute(i, j) = 4;
            y1Face = [y1Face ; P];
        elseif all(Points(P, 3) == 0)
            FaceAttribute(i, j) = 5;
            z0Face = [z0Face ; P];
        elseif all(Points(P, 3) == 1)
            FaceAttribute(i, j) = 6;
            z1Face = [z1Face ; P];
        end
    end
end
NumGlobalFaces = NumGlobalFaces - 1;

figure(1)
subplot(1, 2, 2)
patch('Vertices', Points, 'Faces', x0Face, ...
         'FaceVertexCData', zeros(size(x0Face, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'r'); hold on
patch('Vertices', Points, 'Faces', x1Face, ...
         'FaceVertexCData', zeros(size(x1Face, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'b'); hold on
patch('Vertices', Points, 'Faces', y0Face, ...
         'FaceVertexCData', zeros(size(y0Face, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'g'); hold on
patch('Vertices', Points, 'Faces', y1Face, ...
         'FaceVertexCData', zeros(size(y1Face, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'c'); hold on
patch('Vertices', Points, 'Faces', z0Face, ...
         'FaceVertexCData', zeros(size(z0Face, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'y'); hold on
patch('Vertices', Points, 'Faces', z1Face, ...
         'FaceVertexCData', zeros(size(z1Face, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'magenta'); hold on

A = sparse(NumEles * 4, NumEles * 4);
B = sparse(NumEles * 4, NumEles);
C = sparse(NumEles * 4, NumGlobalFaces);

delta_lm = eye(4);

Volume_eachTet = zeros(NumEles, 1);

for i = 1:NumEles
    V = tetrahedronVolume(Points(Element(i, 1), :), Points(Element(i, 2), :), Points(Element(i, 3), :), Points(Element(i, 4), :));
    Volume_eachTet(i) = V;

    A_loc = zeros(4, 4);
    for j = 1:4
        P1_ID_j = Element(i, FaceLocalOrder(j, 1));
        P2_ID_j = Element(i, FaceLocalOrder(j, 2));
        P3_ID_j = Element(i, FaceLocalOrder(j, 3));

        A_j = Area_tri(Points(P1_ID_j, :), Points(P2_ID_j, :), Points(P3_ID_j, :));
        
        for k = 1:4
            P1_ID_k = Element(i, FaceLocalOrder(k, 1));
            P2_ID_k = Element(i, FaceLocalOrder(k, 2));
            P3_ID_k = Element(i, FaceLocalOrder(k, 3));
            A_k = Area_tri(Points(P1_ID_k, :), Points(P2_ID_k, :), Points(P3_ID_k, :));
            for l = 1:4
                for m = 1:4
                    lambda_lm = V / (60 .* (1 + 2 * delta_lm(l, m)));
                    A_loc(j, k) = A_loc(j, k) + A_j * A_k / (9 * V .^ 2) * dot((Points(Element(i, l), :) - Points(Element(i, j), :)), (Points(Element(i, m), :) - Points(Element(i, k), :))) * lambda_lm;
                end
            end
        end
        B([(i-1) * 4 + j], i) = -A_j;
        C([(i-1) * 4 + j], FaceGlobalID(i, j)) = A_j;
    end
    A([(i-1) * 4 + 1: i * 4], [(i-1) * 4 + 1: i * 4]) = A_loc;
end

K = [A, B, C;
    B', sparse(NumEles, NumEles), sparse(NumEles, NumGlobalFaces);
    C', sparse(NumGlobalFaces, NumEles), sparse(NumGlobalFaces, NumGlobalFaces)];

b = sparse(NumEles * 5 + NumGlobalFaces, 1);

rowsToDelete = [];
for i = 1:NumEles
    for j = 1:4
        FaceID_G = FaceGlobalID(i, j);
        if (FaceAttribute(i, j) == 1 || ...
                FaceAttribute(i, j) == 2 || ...
                FaceAttribute(i, j) == 3 || ...
                FaceAttribute(i, j) == 4)
            if (K(NumEles * 5 + FaceID_G, (i-1) * 4 + j) ~= 0)
                b(NumEles * 5 + FaceID_G) = K(NumEles * 5 + FaceID_G, (i-1) * 4 + j) .* 0; 
                continue;
            else
                fprintf('Wrong 2')
                exit
            end
        end
        
        P1_ID_j = Element(i, FaceLocalOrder(j, 1));
        P2_ID_j = Element(i, FaceLocalOrder(j, 2));
        P3_ID_j = Element(i, FaceLocalOrder(j, 3));

        A_j = Area_tri(Points(P1_ID_j, :), Points(P2_ID_j, :), Points(P3_ID_j, :));

        if (FaceAttribute(i, j) == 5)
            b((i-1) * 4 + j) = b((i-1) * 4 + j) - 1 * A_j; 
            rowsToDelete = [rowsToDelete, FaceID_G];
        elseif (FaceAttribute(i, j) == 6)
            b((i-1) * 4 + j) = b((i-1) * 4 + j) - 100 * A_j; 
            rowsToDelete = [rowsToDelete, FaceID_G];
        end
    end
end

rowsToKeep = setdiff(1:size(K, 1), rowsToDelete + NumEles * 5);
K = K(rowsToKeep, :);
K = K(:, rowsToKeep);
b = b(rowsToKeep);
x = full(inv(K) * b);
pressureEle = x(NumEles * 4 + 1:NumEles * 5);
q_n_face = x(1:NumEles * 4);

figure(2)
subplot(1, 2, 1)
view(3)
title('Show face normal velocity')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', Points, 'Faces', [Element(:, [2, 3, 4]); Element(:, [3, 4, 1]); ...
            Element(:, [4, 1, 2]); ...
            Element(:, [1, 2, 3])], ...
         'FaceVertexCData', zeros(size(Element, 1) * 4, 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0);  hold on
pbaspect([1, 1, 1])

NormalVector_each_face = zeros(4 * NumEles, 3);
Center_eachFace = zeros(4 * NumEles, 3);

Center_eachTet = zeros(NumEles, 3);
q_eachTet = zeros(NumEles, 3);

for i = 1:NumEles
    CentersTET = 1/4.0 * (Points(Element(i, 1), :) + Points(Element(i, 2), :) + Points(Element(i, 3), :) + Points(Element(i, 4), :));
    Center_eachTet(i, :) = CentersTET;

    for j = 1:4
        P1_ID_j = Element(i, FaceLocalOrder(j, 1));
        P2_ID_j = Element(i, FaceLocalOrder(j, 2));
        P3_ID_j = Element(i, FaceLocalOrder(j, 3));

        CenterFace_j = 1.0 / 3.0 .* (Points(P1_ID_j, :) + Points(P2_ID_j, :) + Points(P3_ID_j, :));
        Center_eachFace((i - 1) * 4 + j, :) = CenterFace_j;

        A_j = Area_tri(Points(P1_ID_j, :), Points(P2_ID_j, :), Points(P3_ID_j, :));

        AB = Points(P2_ID_j, :) - Points(P1_ID_j, :);
        AC = Points(P3_ID_j, :) - Points(P1_ID_j, :);
        NormalVector_each_face((i - 1) * 4 + j, :) = cross(AB, AC);
        NormalVector_each_face((i - 1) * 4 + j, :) = NormalVector_each_face((i - 1) * 4 + j, :) ./ norm(NormalVector_each_face((i - 1) * 4 + j, :));

        KO = CentersTET - CenterFace_j;
        KO = KO / norm(KO);

        if (dot(KO, NormalVector_each_face((i - 1) * 4 + j, :) ) > 0)
            NormalVector_each_face((i - 1) * 4 + j, :) = -NormalVector_each_face((i - 1) * 4 + j, :);
        end

        q_eachTet(i, :) = q_eachTet(i, :) + A_j / (3 * Volume_eachTet(i)) .* ( CentersTET - Points(Element(i, j), :)) .* q_n_face((i - 1) * 4 + j);
    end
end

q_n_vector = q_n_face .* NormalVector_each_face;
quiver3(Center_eachFace(:, 1), Center_eachFace(:, 2), Center_eachFace(:, 3), ...
    q_n_vector(:, 1), q_n_vector(:, 2), q_n_vector(:, 3), 'color', 'r');



figure(2)
subplot(1, 2, 2)
view(3)
title('Show tetrahedron center velocity')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', Points, 'Faces', [Element(:, [2, 3, 4]); Element(:, [3, 4, 1]); ...
            Element(:, [4, 1, 2]); ...
            Element(:, [1, 2, 3])], ...
         'FaceVertexCData', zeros(size(Element, 1) * 4, 1), 'FaceColor', 'flat', 'EdgeAlpha', 0.2, 'facealpha', 0);  hold on
pbaspect([1, 1, 1])
quiver3(Center_eachTet(:, 1), Center_eachTet(:, 2), Center_eachTet(:, 3), ...
    q_eachTet(:, 1), q_eachTet(:, 2), q_eachTet(:, 3), 'color', 'r');

figure(3)
view(3)
title('Show pressure')
xlabel('x')
ylabel('y')
zlabel('z')
hold on

tetramesh(Element, Points, pressureEle); hold on; colorbar
pbaspect([1, 1, 1])