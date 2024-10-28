clear all
close all
clc

DFM_mesh;

triangles = msh.TRIANGLES;
points = msh.POS;

tetrahedrons = msh.TETS;
NumTets = size(tetrahedrons, 1);

numTriangles = size(triangles, 1);

figure(1)
view(3)
title('triangles')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', triangles(1:numTriangles, 1:3), 'FaceVertexCData', zeros(numTriangles, 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
pbaspect([1, 1, 1])
% identify frac mesh
% 

Tri_Frac = zeros(numTriangles, 3);
Tri_Bound_N = zeros(numTriangles, 3);
Tri_Bound_D = zeros(numTriangles, 3);

BoundValueTri_N = zeros(numTriangles, 1);
BoundValueTri_D = zeros(numTriangles, 1);

NumFracTri = 1;
NumBoundTri_D = 1;
NumBoundTri_N = 1;

string_N = "Tri_Bound_N(NumBoundTri_N, :) = [ID1, ID2, ID3]; BoundValueTri_N(NumBoundTri_N) = 0; NumBoundTri_N = NumBoundTri_N + 1; ";
string_D_1 = "Tri_Bound_D(NumBoundTri_D, :) = [ID1, ID2, ID3]; BoundValueTri_D(NumBoundTri_D) = 1; NumBoundTri_D = NumBoundTri_D + 1;";
string_D_0 = "Tri_Bound_D(NumBoundTri_D, :) = [ID1, ID2, ID3]; BoundValueTri_D(NumBoundTri_D) = 0; NumBoundTri_D = NumBoundTri_D + 1;";

for i = 1:numTriangles
    ID1 = triangles(i, 1);
    ID2 = triangles(i, 2);
    ID3 = triangles(i, 3);

    if (points(ID1, 1) == 0 && ...
        points(ID2, 1) == 0 && ...
        points(ID3, 1) == 0)
        eval(string_N);
    elseif (points(ID1, 1) == 1 && ...
        points(ID2, 1) == 1 && ...
        points(ID3, 1) == 1)
        eval(string_N);
    elseif (points(ID1, 2) == 0 && ...
            points(ID2, 2) == 0 && ...
            points(ID3, 2) == 0)
        eval(string_N);
    elseif (points(ID1, 2) == 1 && ...
            points(ID2, 2) == 1 && ...
            points(ID3, 2) == 1)
        eval(string_N);
    elseif (points(ID1, 3) == 0 && ...
            points(ID2, 3) == 0 && ...
            points(ID3, 3) == 0)
        eval(string_D_0);
    elseif (points(ID1, 3) == 1 && ...
            points(ID2, 3) == 1 && ...
            points(ID3, 3) == 1)
        eval(string_D_1);

    else
        
        Tri_Frac(NumFracTri, :) = [ID1, ID2, ID3];
        NumFracTri = NumFracTri + 1;
    end
end

NumFracTri = NumFracTri - 1;
NumBoundTri_N = NumBoundTri_N - 1;
NumBoundTri_D = NumBoundTri_D - 1;

Tri_Frac(NumFracTri + 1:end, :) = [];
Tri_Bound_N(NumBoundTri_N + 1:end, :) = [];
Tri_Bound_D(NumBoundTri_D + 1:end, :) = [];

BoundValueTri_N(NumBoundTri_N + 1:end, :) = [];
BoundValueTri_D(NumBoundTri_D + 1:end, :) = [];

figure(2)
view(3)
title('diff triangles')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', Tri_Frac, 'FaceVertexCData', zeros(NumFracTri, 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'r'); hold on
patch('Vertices', points, 'Faces', Tri_Bound_N, 'FaceVertexCData', zeros(NumBoundTri_N, 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'b'); hold on
patch('Vertices', points, 'Faces', Tri_Bound_D, 'FaceVertexCData', zeros(NumBoundTri_D, 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'g'); hold on
pbaspect([1, 1, 1])

Tri_Frac_sort = sort(Tri_Frac, 2);
Tri_Bound_N_sort = sort(Tri_Bound_N, 2);
Tri_Bound_D_sort = sort(Tri_Bound_D, 2);

%---------------------------------

Info_tet = zeros(NumTets * 4, 4);
% 1: global numbering
% 2: attri 0-common_interface, 1-frac, 2-D, 3-N
% 3: Frac ele No if it is
% 4: sign
% note: when it is a frac face, it is isolated to a discontinous one, having a unique
% numbering ID
% The tet are with MHFEM, only frac triangles are hybridized

Info_tri = zeros(NumFracTri * 3, 3);
% 1: global edge numbering ID
% 2: attri 0-common_interface, 1-nothing, 2-D, 3-N
% 3: D or N boundary condition value if it is D or N
% note: frac triangles are all with MHFEM, all interior edges are
% hybridized


FaceLocalOrder = [2, 3, 4
   3, 4, 1
   4, 1, 2
   1, 2, 3];
EdgeLocalOrder = [2, 3
    3, 1
    1, 2];

%--------------edges
%--------------edges
%--------------edges
% temp variable
EdgeAcc = zeros(NumFracTri * 3, 3);

NumInteriosEdge = 1;
NumNeumannEdge = 1;
NumDirichiletEdge = 1;

NumGlobalEdges = 1;
for i = 1:NumFracTri
    for j = 1:3
        edgeLocalNo = (i-1) * 3 + j;
        ID1 = Tri_Frac(i, EdgeLocalOrder(j, 1)); 
        ID2 = Tri_Frac(i, EdgeLocalOrder(j, 2)); 
        
        ID_sort = sort([ID1, ID2]);

        [isRowPresent, rowIndex] = ismember(ID_sort, EdgeAcc(1:NumGlobalEdges, 1:2), 'rows');

        if ~isRowPresent
            Info_tri(edgeLocalNo, 1) = NumGlobalEdges;

            EdgeAcc(NumGlobalEdges, :) = [ID_sort, edgeLocalNo];
            NumGlobalEdges = NumGlobalEdges + 1;
        else
            edgeLocalNo_II = EdgeAcc(rowIndex, 3);

            Info_tri(edgeLocalNo, :) = Info_tri(edgeLocalNo_II, :);
            continue
        end

        % boundary identification
        % boundary identification
        % boundary identification triangles
        if (points(ID1, 3) == 0 && points(ID2, 3) == 0)
            Info_tri(edgeLocalNo, 2) = 2;
            Info_tri(edgeLocalNo, 3) = 0;
            NumDirichiletEdge = NumDirichiletEdge + 1;
            continue;
        end

        if (points(ID1, 3) == 1 && points(ID2, 3) == 1)
            Info_tri(edgeLocalNo, 2) = 2;
            Info_tri(edgeLocalNo, 3) = 1;
            NumDirichiletEdge = NumDirichiletEdge + 1;
            continue;
        end
        
        if (points(ID1, 1) == 0 && points(ID2, 1) == 0) || ...
           (points(ID1, 1) == 1 && points(ID2, 1) == 1) || ...
           (points(ID1, 2) == 0 && points(ID2, 2) == 0) || ...
           (points(ID1, 2) == 1 && points(ID2, 2) == 1)
            Info_tri(edgeLocalNo, 2) = 3;
            Info_tri(edgeLocalNo, 3) = 0.;
            NumNeumannEdge = NumNeumannEdge + 1;
            continue;
        end
        
        % common interface edge
        Info_tri(edgeLocalNo, 2) = 0;
        Info_tri(edgeLocalNo, 3) = 0.;
        NumInteriosEdge = NumInteriosEdge + 1;

    end
end
clear EdgeAcc
NumGlobalEdges = NumGlobalEdges - 1;
NumInteriosEdge = NumInteriosEdge - 1;
NumNeumannEdge = NumNeumannEdge - 1;
NumDirichiletEdge = NumDirichiletEdge - 1;

%------------triangles
%------------triangles
%------------triangles

% temp variable
TriAcc = zeros(NumTets * 4, 4);

NumGlobalTri = 1;
for i = 1:NumTets
    for j = 1:4
        triangleLocalNo = (i-1) * 4 + j;

        ID1 = tetrahedrons(i, FaceLocalOrder(j, 1)); 
        ID2 = tetrahedrons(i, FaceLocalOrder(j, 2)); 
        ID3 = tetrahedrons(i, FaceLocalOrder(j, 3));
        
        ID_sort = sort([ID1, ID2, ID3]);
        
        [isRowPresent, rowIndex] = ismember(ID_sort, TriAcc(1:NumGlobalTri, 1:3), 'rows');

        if ~isRowPresent
            Info_tet(triangleLocalNo, 1) = NumGlobalTri;
            Info_tet(triangleLocalNo, 4) = 1;

            TriAcc(NumGlobalTri, :) = [ID_sort, triangleLocalNo];

            NumGlobalTri = NumGlobalTri + 1;
        else
            triangleLocalNo_II = TriAcc(rowIndex, 4);

            Info_tet(triangleLocalNo, :) = Info_tet(triangleLocalNo_II, :);
            Info_tet(triangleLocalNo, 4) = -Info_tet(triangleLocalNo, 4);
            if (Info_tet(triangleLocalNo, 2) == 1) % it is a frac triangle
                Info_tet(triangleLocalNo, 1) = NumGlobalTri;
                Info_tet(triangleLocalNo, 4) = 1;
                NumGlobalTri = NumGlobalTri + 1;
            end
            continue
        end

        [isRowPresent, rowIndex] = ismember(ID_sort, Tri_Frac_sort, 'rows');

        if isRowPresent
            Info_tet(triangleLocalNo, 2) = 1;
            Info_tet(triangleLocalNo, 3) = rowIndex; % the ID of the frac element No
            continue
        end

        [isRowPresent, rowIndex] = ismember(ID_sort, Tri_Bound_D_sort, 'rows');
        if isRowPresent
            Info_tet(triangleLocalNo, 2) = 2;
            continue
        end

        [isRowPresent, rowIndex] = ismember(ID_sort, Tri_Bound_N_sort, 'rows');
        if isRowPresent
            Info_tet(triangleLocalNo, 2) = 3;
            continue
        end
        Info_tet(triangleLocalNo, 2) = 0;
    end
end
NumGlobalTri = NumGlobalTri - 1;
clear TriAcc

%--------------check the numbering system
%--------------check the numbering system
%--------------check the numbering system
tets_frac_adjacent = find(Info_tet(:, 2) == 1);
tets_frac_adjacent = ceil(tets_frac_adjacent./4);

tets_D_adjacent = find(Info_tet(:, 2) == 2);
tets_D_adjacent = ceil(tets_D_adjacent./4);

tets_N_adjacent = find(Info_tet(:, 2) == 3);
tets_N_adjacent = ceil(tets_N_adjacent./4);

figure(3)
subplot(1, 3, 1)
view(3)
title('frac releted tets')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetramesh(tetrahedrons(tets_frac_adjacent, [1:4]), points, tetrahedrons(tets_frac_adjacent, 1) .* 0, 'FaceAlpha', 1., ...
    'Edgecolor', 'g'); hold on;
pbaspect([1, 1, 1])

subplot(1, 3, 2)
view(3)
title('D releted tets')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetramesh(tetrahedrons(tets_D_adjacent, [1:4]), points, tetrahedrons(tets_D_adjacent, 1) .* 0, 'FaceAlpha', 1., ...
    'Edgecolor', 'b'); hold on;
pbaspect([1, 1, 1])

subplot(1, 3, 3)
view(3)
title('N releted tets')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetramesh(tetrahedrons(tets_N_adjacent, [1:4]), points, tetrahedrons(tets_N_adjacent, 1) .* 0, 'FaceAlpha', 1., ...
    'Edgecolor', 'r'); hold on;
pbaspect([1, 1, 1])

%--------edges checking numbering system
%--------edges checking numbering system
%--------edges checking numbering system

tets_D_adjacent = find(Info_tri(:, 2) == 2);
tets_D_adjacent = ceil(tets_D_adjacent./3);

tets_N_adjacent = find(Info_tri(:, 2) == 3);
tets_N_adjacent = ceil(tets_N_adjacent./3);

tets_I_adjacent = find(Info_tri(:, 2) == 0);
tets_I_adjacent = ceil(tets_I_adjacent./3);
tets_I_adjacent = unique(tets_I_adjacent);

figure(4)
subplot(1, 3, 1)
view(3)
title('interior tri')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', Tri_Frac(tets_I_adjacent, 1:3), 'FaceVertexCData', zeros(size(tets_I_adjacent, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
pbaspect([1, 1, 1])

subplot(1, 3, 2)
view(3)
title('d tri')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', Tri_Frac(tets_D_adjacent, 1:3), 'FaceVertexCData', zeros(size(tets_D_adjacent, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
pbaspect([1, 1, 1])

subplot(1, 3, 3)
view(3)
title('n tri')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', Tri_Frac(tets_N_adjacent, 1:3), 'FaceVertexCData', zeros(size(tets_N_adjacent, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
pbaspect([1, 1, 1])

clear tets_frac_adjacent tets_D_adjacent tets_N_adjacent