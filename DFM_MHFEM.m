%clear all
close all
clc

% please run AddressingDFMBoundary first
if ~exist('Info_tet', 'var')
    clear all
    disp('run AddressingDFMBoundary ...')
    AddressingDFMBoundary;
else
    disp('AddressingDFMBoundary has been run...')
end
close all

Permeability_Tensor = eye(3) .* 1e-5;
Permeability_Tensor_inv = inv(Permeability_Tensor);

Frac_conductivity = 1e-3;

Dim = NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + ...
    NumGlobalEdges;

K = sparse(Dim, Dim);
b_right = sparse(Dim, 1);

Volume_eachTet = zeros(NumTets, 1);
delta_lm = eye(4);
%sign_check = zeros(NumGlobalTri, 2);

for ele = 1:NumTets
    V = tetrahedronVolume(points(tetrahedrons(ele, 1), :), ...
        points(tetrahedrons(ele, 2), :), points(tetrahedrons(ele, 3), :), points(tetrahedrons(ele, 4), :));
    Volume_eachTet(ele) = V;
    
    A_loc = zeros(4, 4);
    B_loc = zeros(4, 1);
    globalTriID = zeros(4, 1);

    for i = 1:4
        NodeFace_i = [tetrahedrons(ele, FaceLocalOrder(i, 1));
                tetrahedrons(ele, FaceLocalOrder(i, 2));
                tetrahedrons(ele, FaceLocalOrder(i, 3));];
        A_i = Area_tri(points(NodeFace_i(1), :), points(NodeFace_i(2), :), points(NodeFace_i(3), :));
        
        localID_i = (ele - 1) * 4 + i;
        globalTriID(i) = Info_tet(localID_i, 1);
        sign_i = Info_tet(localID_i, 4);
        
        % for p = 1:2
        %     if sign_check(globalTriID(i), p) == 0
        %         sign_check(globalTriID(i), p) = sign_i;
        %         break
        %     end
        %     if p == 2
        %         error('wrong')
        %     end
        % end

        for j = 1:4
            NodeFace_j = [tetrahedrons(ele, FaceLocalOrder(j, 1));
                tetrahedrons(ele, FaceLocalOrder(j, 2));
                tetrahedrons(ele, FaceLocalOrder(j, 3));];
            A_j = Area_tri(points(NodeFace_j(1), :), points(NodeFace_j(2), :), points(NodeFace_j(3), :));

            localID_j = (ele - 1) * 4 + j;
            sign_j = Info_tet(localID_j, 4);

            for k = 1:3       
                for l = 1:3  
                    for a = 1:4
                        for b = 1:4
                            A_loc(i, j) = A_loc(i, j) + sign_i * sign_j * A_i * A_j / (9 * V^2) * ...
                                Permeability_Tensor_inv(k, l) * (points(tetrahedrons(ele, a), k) - points(tetrahedrons(ele, i), k)) * (points(tetrahedrons(ele, b), l) - ...
                                points(tetrahedrons(ele, j), l)) * V / 20 .* (1 + 1 * delta_lm(a, b));
                        end
                        
                    end
                end
            end
        end
        B_loc(i) = sign_i * A_i;

        if Info_tet(localID_i, 2) == 1
            corresponding_frac_tri_ID = Info_tet(localID_i, 3);
            %disp(['corresponding_frac_tri_ID = ', num2str(corresponding_frac_tri_ID)]);

            K(globalTriID(i), NumGlobalTri + NumTets + NumFracTri * 3 + corresponding_frac_tri_ID) = A_i;
            K(NumGlobalTri + NumTets + NumFracTri * 3 + corresponding_frac_tri_ID, globalTriID(i)) = A_i;
        end

        if (Info_tet(localID_i, 2) == 2) % Dirichilet boundary
            b_right(globalTriID(i), 1) = b_right(globalTriID(i), 1) - sign_i * Info_tet(localID_i, 3) * A_i;
        end

    end
    
    K(globalTriID, globalTriID) = K(globalTriID, globalTriID) + A_loc;

    K(globalTriID, ele + NumGlobalTri) = K(globalTriID, ele + NumGlobalTri) - B_loc;

    K(ele + NumGlobalTri, globalTriID) = K(ele + NumGlobalTri, globalTriID) - B_loc';
end

for i = 1:size(Info_tet, 1)
    if (Info_tet(i, 2) == 3)
        % Neumann BC
        globalTriID = Info_tet(i, 1);
        b_right = b_right + K(:, globalTriID) * Info_tet(i, 3);
        %disp(nnz(K(globalTriID, :)))
        K(globalTriID, :) = 0;
        K(:, globalTriID) = 0;
        K(globalTriID,globalTriID) = 1;
        b_right(globalTriID) = Info_tet(i, 3) * 1;
    end
    
end

%-----------------------------frac triangles MHFEM-----------
%-----------------------------frac triangles MHFEM-----------
%-----------------------------frac triangles MHFEM-----------
M =     [[2., 0., 1., 0., 1., 0.];
    [0., 2., 0., 1., 0., 1.];
    [1., 0., 2., 0., 1., 0.];
    [0., 1., 0., 2., 0., 1.];
    [1., 0., 1., 0., 2., 0.];
    [0., 1., 0., 1., 0., 2.]];
for ele = 1:NumFracTri
    P1 = points(Tri_Frac(ele, 1), :);
    P2 = points(Tri_Frac(ele, 2), :);
    P3 = points(Tri_Frac(ele, 3), :);

    if P1(1) == 0.5 && P2(1) == 0.5 && P3(1) == 0.5
        P1(1) = []; P2(1) = []; P3(1) = [];
    elseif P1(3) == 0.5 && P2(3) == 0.5 && P3(3) == 0.5
        P1(3) = []; P2(3) = []; P3(3) = [];
    else
        error("errous frac tri")
    end

    T_area = Area_tri(P1, P2, P3);
    
    N = zeros(6, 3);
    N(2 + 1 :4, 1 + 0) = (P2 - P1)';
    N(4 + 1 :6, 1 + 0) = (P3 - P1)';
    N(0 + 1 :2, 1 + 1) = (P1 - P2)';
    N(4 + 1 :6, 1 + 1) = (P3 - P2)';
    N(0 + 1 :2, 1 + 2) = (P1 - P3)';
    N(2 + 1 :4, 1 + 2) = (P2 - P3)';

    C = [
        [norm(P3-P2), 0, 0],
        [0, norm(P3-P1), 0],
        [0, 0, norm(P1-P2)]
    ];

    A_loc = 1 ./ Frac_conductivity .* 1. / 48 / T_area * C' * N' * M * N * C;
    
    K(((ele-1)*3+1:ele*3) + NumGlobalTri + NumTets, ((ele-1)*3+1:ele*3) + NumGlobalTri + NumTets) = A_loc;

    K(((ele-1)*3+1:ele*3) + NumGlobalTri + NumTets, NumGlobalTri + NumTets + NumFracTri * 3 + ele) = -diag(C);

    K(NumGlobalTri + NumTets + NumFracTri * 3 + ele, ((ele-1)*3+1:ele*3) + NumGlobalTri + NumTets) = -diag(C)';
    

    % globalEdgeID_tt = Info_tri((ele-1) * 3 + 1 : (ele-1) * 3 + 3, 1);

    for i = 1:3
        globalEdgeID_tt = Info_tri((ele-1) * 3 + i, 1);
        K(((ele-1)*3+i) + NumGlobalTri + NumTets, ...
            NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt) ...
            = C(i, i);
        K(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt, ...
            ((ele-1)*3+i) + NumGlobalTri + NumTets) ...
            = C(i, i);

        if (Info_tri((ele-1) * 3 + i, 2) == 2)
            b_right(((ele-1)*3+i) + NumGlobalTri + NumTets) = -Info_tri((ele-1) * 3 + i, 3) * C(i, i);
        end
    end

end

for i = 1:NumFracTri*3
    if (Info_tri(i, 2) == 3) % Neumann
        globalEdgeID_tt = Info_tri(i, 1);

        %b_right = b_right-K(:, NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt) * Info_tri(i, 3);
        b_right(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt) =...
            K(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt, ...
            NumGlobalTri + NumTets  + i) * Info_tri(i, 3);
        %b_right(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt)
        % K(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt, :) = 0;
        % K(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt, NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt) = 1;
        % b_right(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt) = Info_tri(i, 3);
    end
end
for i = 1:NumFracTri*3
    if (Info_tri(i, 2) == 2)
        globalEdgeID_tt = Info_tri(i, 1);
        K(:, NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt) = 0;
        K(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt, :) = 0;
        b_right(NumGlobalTri + NumTets + NumFracTri * 3 + NumFracTri + globalEdgeID_tt) = 0;
    end
end
% -----------------------------------------------------------

x_fc = full((K(NumGlobalTri + NumTets + 1:end, NumGlobalTri + NumTets + 1:end)) \ b_right(NumGlobalTri + NumTets + 1:end));
pressure_pure_fc = x_fc(NumFracTri * 3 + 1: NumFracTri * 4);

figure(4);
subplot(1, 3, 1)
spy(K, 'k.'); 
title('Sparsity Pattern');
xlabel('Column Index');
ylabel('Row Index');

subplot(1, 3, 2)
title('Sparsity Pattern');
spy([K(NumGlobalTri + NumTets + 1:end, NumGlobalTri + NumTets + 1:end), b_right(NumGlobalTri + NumTets + 1:end)], 'k.'); 
xlabel('Column Index');
ylabel('Row Index');

subplot(1, 3, 3)
view(3)
title('Pressure pure frac')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', Tri_Frac, 'FaceVertexCData', pressure_pure_fc, 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 1); hold on
colorbar
pbaspect([1, 1, 1])


x_couple = full((K) \ b_right);
pressure_tet = x_couple(NumGlobalTri + 1:NumGlobalTri + NumTets);
pressure_tri = x_couple(NumGlobalTri + NumTets + NumFracTri * 3 + 1:NumGlobalTri + NumTets + NumFracTri * 4);
figure(5)
subplot(1, 2, 1)
view(3)
title('Show pressure tet coupled')
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetramesh(tetrahedrons, points, pressure_tet); hold on; colorbar
pbaspect([1, 1, 1])

subplot(1, 2, 2)
view(3)
title('Show pressure tri coupled')
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', Tri_Frac, ...
    'FaceVertexCData', pressure_tri, 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 1); hold on
colorbar
pbaspect([1, 1, 1])

return 


pressureEle = x(NumGlobalTri + 1:NumGlobalTri + NumTets);

% figure(6)
% view(3)
% title('Check numbering system')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% 
% tetramesh(tetrahedrons, points, pressureEle, 'facealpha', 0); hold on
% 
% visited=zeros(NumGlobalTri, 1);
% for ele = 1:NumTets
%     for i = 1:4
%         NodeFace_i = [tetrahedrons(ele, FaceLocalOrder(i, 1));
%                 tetrahedrons(ele, FaceLocalOrder(i, 2));
%                 tetrahedrons(ele, FaceLocalOrder(i, 3));];
%         P_center = 1/3. .* (points(NodeFace_i(1), :) + ... 
%             points(NodeFace_i(2), :) + ...
%             points(NodeFace_i(3), :));
% 
%         localID_i = (ele - 1) * 4 + i;
% 
%         if visited(Info_tet(localID_i, 1)) == 0
%             text(P_center(1), P_center(2), P_center(3), num2str(Info_tet(localID_i, 1)), 'Color', 'k');
%         else
%             text(P_center(1), P_center(2), P_center(3), num2str(Info_tet(localID_i, 1)), 'Color', 'r');
%         end
% 
%         visited(Info_tet(localID_i, 1)) = 1;
%     end
% end
% 
