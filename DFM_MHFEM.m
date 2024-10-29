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

Permeability_Tensor = eye(3);
Permeability_Tensor_inv = inv(Permeability_Tensor);

Dim = NumGlobalTri + NumTets + NumFracTri;%NumFracTri * 3 + NumFracTri + ...
    %NumInteriosEdge;

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

            K(globalTriID(i), NumGlobalTri + NumTets + corresponding_frac_tri_ID) = A_i;
            K(NumGlobalTri + NumTets + corresponding_frac_tri_ID, globalTriID(i)) = A_i;
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
        % nnz(K(globalTriID, :))
        K(globalTriID, :) = 0;
        K(:, globalTriID) = 0;
        K(globalTriID,globalTriID) = 1;
        b_right(globalTriID) = Info_tet(i, 3) * 1;
    end
    
end

figure(4);
spy(K, 'k.'); 
title('Sparsity Pattern');
xlabel('Column Index');
ylabel('Row Index');

x = full((K) \ b_right);
pressureEle = x(NumGlobalTri + 1:NumGlobalTri + NumTets);
figure(5)
view(3)
title('Show pressure')
xlabel('x')
ylabel('y')
zlabel('z')
hold on

tetramesh(tetrahedrons, points, pressureEle); hold on; colorbar

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
