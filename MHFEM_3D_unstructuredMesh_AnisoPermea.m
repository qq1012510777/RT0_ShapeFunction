clc
clear all
close all

box_tets_3D;

% Permeability_Tensor = [1, 3, 1;
%     1, 2, 4;
%     0, 0, 3];
% Permeability_Tensor = (Permeability_Tensor + Permeability_Tensor') / 2

% for i = 1:100
%     Permeability_Tensor = rand(3, 3);
%     Permeability_Tensor = (Permeability_Tensor + Permeability_Tensor') / 2;
%     Permeability_Tensor = Permeability_Tensor + 3 * eye(3);
% 
%     is_pos_def = all(eig(Permeability_Tensor) > 0);
%     if (is_pos_def == 0)
%         continue
%     else
%         fprintf("Permeability tensor is:\n")
%         disp(Permeability_Tensor)
%         break;
%     end
% end

Permeability_Tensor = [    3.7922    0.4976    0.6672
    0.4976    3.8491    0.8459
    0.6672    0.8459    3.7431]

% Permeability_Tensor = [    3.7922    0 0
%     0   3.8491    0
%     0   0   3.7431]

Permeability_Tensor_inv = inv(Permeability_Tensor);

Points = msh.POS;
Element = msh.TETS(:, 1:4);

min_x = min(Points(:, 1));
max_x = max(Points(:, 1));
min_y = min(Points(:, 2));
max_y = max(Points(:, 2));
min_z = min(Points(:, 3));
max_z = max(Points(:, 3));

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
        
        if all(Points(P, 1) == min_x)
            FaceAttribute(i, j) = 1;
            x0Face = [x0Face ; P];
        elseif all(Points(P, 1) == max_x)
            FaceAttribute(i, j) = 2;
            x1Face = [x1Face ; P];
        elseif all(Points(P, 2) == min_y)
            FaceAttribute(i, j) = 3;
            y0Face = [y0Face ; P];
        elseif all(Points(P, 2) == max_y)
            FaceAttribute(i, j) = 4;
            y1Face = [y1Face ; P];
        elseif all(Points(P, 3) == min_z)
            FaceAttribute(i, j) = 5;
            z0Face = [z0Face ; P];
        elseif all(Points(P, 3) == max_z)
            FaceAttribute(i, j) = 6;
            z1Face = [z1Face ; P];
        end
    end
end
NumGlobalFaces = NumGlobalFaces - 1;

A = sparse(NumEles * 4, NumEles * 4);
B = sparse(NumEles * 4, NumEles);
C = sparse(NumEles * 4, NumGlobalFaces);

delta_lm = eye(4);
Volume_eachTet = zeros(NumEles, 1);

Q_sss = zeros(3, 3);
for dir = 1:3
    for ele = 1:NumEles
        V = tetrahedronVolume(Points(Element(ele, 1), :), Points(Element(ele, 2), :), Points(Element(ele, 3), :), Points(Element(ele, 4), :));
        Volume_eachTet(ele) = V;
    
        A_loc = zeros(4, 4);
        for i = 1:4
            NodeFace_i = [Element(ele, FaceLocalOrder(i, 1));
                          Element(ele, FaceLocalOrder(i, 2));
                          Element(ele, FaceLocalOrder(i, 3));];
            A_i = Area_tri(Points(NodeFace_i(1), :), Points(NodeFace_i(2), :), Points(NodeFace_i(3), :));
        
            for j = 1:4
                NodeFace_j = [Element(ele, FaceLocalOrder(j, 1));
                              Element(ele, FaceLocalOrder(j, 2));
                              Element(ele, FaceLocalOrder(j, 3));];
                A_j = Area_tri(Points(NodeFace_j(1), :), Points(NodeFace_j(2), :), Points(NodeFace_j(3), :));
        
                for k = 1:3       
                    for l = 1:3  
        
                        for a = 1:4
                            for b = 1:4
                                A_loc(i, j) = A_loc(i, j) + A_i * A_j / (9 * V^2) * ...
                                    Permeability_Tensor_inv(k, l) * (Points(Element(ele, a), k) - Points(Element(ele, i), k)) * (Points(Element(ele, b), l) - Points(Element(ele, j), l)) * V / 20 .* (1 + 1 * delta_lm(a, b));
                            end
                            
                        end
                    end
                end
            end
            B([(ele-1) * 4 + i], ele) = -A_i;
            C([(ele-1) * 4 + i], FaceGlobalID(ele, i)) = A_i;
        end
        A([(ele-1) * 4 + 1: ele * 4], [(ele-1) * 4 + 1: ele * 4]) = A_loc;
    end
    
    K = [A, B, C;
        B', sparse(NumEles, NumEles), sparse(NumEles, NumGlobalFaces);
        C', sparse(NumGlobalFaces, NumEles), sparse(NumGlobalFaces, NumGlobalFaces)];
    
    b = sparse(NumEles * 5 + NumGlobalFaces, 1);
    
    rowsToDelete = [];
    if dir == 1
        NeumannBCNo = [3, 4, 5, 6];
        DirichletBCNo = [1, 2];
    elseif dir == 2
        NeumannBCNo = [1, 2, 5, 6];
        DirichletBCNo = [3, 4];
    elseif dir == 3
        NeumannBCNo = [1, 2, 3, 4];
        DirichletBCNo = [5, 6];
    end
    
    for i = 1:NumEles
        for j = 1:4
            FaceID_G = FaceGlobalID(i, j);
            if (FaceAttribute(i, j) == NeumannBCNo(1) || ...
                    FaceAttribute(i, j) == NeumannBCNo(2) || ...
                    FaceAttribute(i, j) == NeumannBCNo(3) || ...
                    FaceAttribute(i, j) == NeumannBCNo(4))
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
    
            if (FaceAttribute(i, j) == DirichletBCNo(1))
                b((i-1) * 4 + j) = b((i-1) * 4 + j) - 0 * A_j; 
                rowsToDelete = [rowsToDelete, FaceID_G];
            elseif (FaceAttribute(i, j) == DirichletBCNo(2))
                b((i-1) * 4 + j) = b((i-1) * 4 + j) - 1 * A_j; 
                rowsToDelete = [rowsToDelete, FaceID_G];
            end
   
    
            
        end
    end
    
    rowsToKeep = setdiff(1:size(K, 1), rowsToDelete + NumEles * 5);
    K = K(rowsToKeep, :);
    K = K(:, rowsToKeep);
    b = b(rowsToKeep);
    x = full((K) \ b);
    pressureEle = x(NumEles * 4 + 1:NumEles * 5);
    q_n_face = x(1:NumEles * 4);
    
    
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

    Q_ppp = q_eachTet .* Volume_eachTet;
    V_dmain = sum(Volume_eachTet);
    Qx = sum(Q_ppp(:, 1)) / V_dmain;
    Qy = sum(Q_ppp(:, 2)) / V_dmain;
    Qz = sum(Q_ppp(:, 3)) / V_dmain;

    Q_sss(:, dir) = [Qx; Qy; Qz];
end
-Q_sss

return 
q_n_vector = q_n_face .* NormalVector_each_face;

figure(2)
view(3)
title('Show tetrahedron center velocity')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetramesh(Element, Points, pressureEle .* 0, 'FaceAlpha', 0.); hold on; 
pbaspect([max_x-min_x, max_y-min_y, max_z-min_z])
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
pbaspect([max_x-min_x, max_y-min_y, max_z-min_z])
    
    % Permeability_Tensor = [    3.7922    0.4976    0.6672
    %                             0.4976    3.8491    0.8459
    %                             0.6672    0.8459    3.7431]

    % l = 0.1 
    % 3.7068    0.2066    0.2920
    % 0.2057    3.7262    0.3792
    % 0.2992    0.3899    3.5934

    % l = 0.05
    % 3.7088    0.2133    0.3007
    % 0.2125    3.7291    0.3897
    % 0.3080    0.4007    3.5968