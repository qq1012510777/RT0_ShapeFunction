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
IfBoundaryTriangle = zeros(numTriangles, 1);

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


%---------------------------------
% 1. global face No, 2. sign for each face, 3. face attribute
% 4. tet no.    5. local face no.       6. FracTri no. or D no.
% or N no. 
GlobalNum_sign_attr = zeros(NumTets * 4, 6);



Tri_Frac_sort = sort(Tri_Frac, 2);
Tri_Bound_N_sort = sort(Tri_Bound_N, 2);
Tri_Bound_D_sort = sort(Tri_Bound_D, 2);

FaceLocalOrder = [2, 3, 4
   3, 4, 1
   4, 1, 2
   1, 2, 3];

globalTriID = 1;
Triang_acc = zeros(NumTets * 4, 4);

for i = 1:NumTets
    disp(['i = ', num2str(i), ', NumTets = ', num2str(NumTets)])
    for j = 1:4
        

        edgeLocalNo = (i-1) * 4 + j;

        ID1 = tetrahedrons(i, FaceLocalOrder(j, 1)); 
        ID2 = tetrahedrons(i, FaceLocalOrder(j, 2)); 
        ID3 = tetrahedrons(i, FaceLocalOrder(j, 3));
        
        ID_sort = sort([ID1, ID2, ID3]);

        [isRowPresent, rowIndex] = ismember(ID_sort, Triang_acc(1:globalTriID, 1:3), 'rows');
        
        GlobalNum_sign_attr(edgeLocalNo, [4, 5]) = [i, j];

        if ~isRowPresent
            GlobalNum_sign_attr(edgeLocalNo, 1) = globalTriID;
            Triang_acc(globalTriID, :) = [ID_sort, edgeLocalNo];
            globalTriID = globalTriID + 1;
        else
            GlobalNum_sign_attr(edgeLocalNo, :) = GlobalNum_sign_attr(Triang_acc(rowIndex, 4), :);
            GlobalNum_sign_attr(edgeLocalNo, 2) = -GlobalNum_sign_attr(edgeLocalNo, 2);
            GlobalNum_sign_attr(edgeLocalNo, [4, 5]) = [i, j];
            continue
        end
        
        [isRowPresent, rowIndex] = ismember(ID_sort, Tri_Frac_sort, 'rows');
        
        if isRowPresent % it is a frac triangle
            GlobalNum_sign_attr(edgeLocalNo, 2) = 1; % sign
            GlobalNum_sign_attr(edgeLocalNo, 3) = 1; 
            GlobalNum_sign_attr(edgeLocalNo, 6) = rowIndex; % tri No.
            continue
        end

        [isRowPresent, rowIndex] = ismember(ID_sort, Tri_Bound_D_sort, 'rows');
        
        if isRowPresent % it is a dirichilet triangle
            GlobalNum_sign_attr(edgeLocalNo, 2) = 1; % sign
            GlobalNum_sign_attr(edgeLocalNo, 3) = 2; 
            GlobalNum_sign_attr(edgeLocalNo, 6) = rowIndex; % D No.
            continue
        end

        [isRowPresent, rowIndex] = ismember(ID_sort, Tri_Bound_N_sort, 'rows');

        if isRowPresent % it is a dirichilet triangle
            GlobalNum_sign_attr(edgeLocalNo, 2) = 1; % sign
            GlobalNum_sign_attr(edgeLocalNo, 3) = 3; 
            GlobalNum_sign_attr(edgeLocalNo, 6) = rowIndex; % N No.
            continue
        end
        
        % then it is just a normal interface triangle
        GlobalNum_sign_attr(edgeLocalNo, 2) = 1; % sign
        GlobalNum_sign_attr(edgeLocalNo, 3) = 0; 
        
    end
end

globalTriID = globalTriID - 1;

tets_Frac_adjacent = find(GlobalNum_sign_attr(:, 3) == 1);
tets_Frac_adjacent = GlobalNum_sign_attr(tets_Frac_adjacent, 4);

figure(3)
view(3)
title('Frac releted tets')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetramesh(tetrahedrons(tets_Frac_adjacent, [1:4]), points, tetrahedrons(tets_Frac_adjacent, 1) .* 0, 'FaceAlpha', 1., ...
    'Edgecolor', 'r'); hold on;
pbaspect([1, 1, 1])

tets_D_adjacent = find(GlobalNum_sign_attr(:, 3) == 2);
tets_D_adjacent = GlobalNum_sign_attr(tets_D_adjacent, 4);

figure(4)
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

tets_N_adjacent = find(GlobalNum_sign_attr(:, 3) == 3);
tets_N_adjacent = GlobalNum_sign_attr(tets_N_adjacent, 4);

figure(5)
view(3)
title('N releted tets')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
tetramesh(tetrahedrons(tets_N_adjacent, [1:4]), points, tetrahedrons(tets_N_adjacent, 1) .* 0, 'FaceAlpha', 1., ...
    'Edgecolor', 'g'); hold on;
pbaspect([1, 1, 1])

%--------------------------------
%--------------------------------
%--------------------------------
NumInteriosEdge = 1;
NumNeumannEdge = 1;
NumDirichiletEdge = 1;

EdgeLocalOrder = [2, 3
    3, 1
    1, 2];

EdgeAcc = zeros(NumFracTri * 3, 5);

EdgeAttr = zeros(NumFracTri * 3, 3);

GlobalEdgeID = 1;

for i = 1:NumFracTri
    for j = 1:3
        edgeLocalID = (i - 1) * 3 + j;

        ID1 = Tri_Frac(i, EdgeLocalOrder(j, 1));
        ID2 = Tri_Frac(i, EdgeLocalOrder(j, 2));

        ID_sort = sort([ID1, ID2]);
        
        [isRowPresent, rowIndex] = ismember(ID_sort, EdgeAcc(1:GlobalEdgeID, 1:2), 'rows');

        if ~isRowPresent
            EdgeAttr(edgeLocalID, 1) = GlobalEdgeID;
      
            EdgeAcc(GlobalEdgeID, :) = [ID_sort, i, j, edgeLocalID];
            GlobalEdgeID = GlobalEdgeID + 1;
        else
            
            EdgeAttr(edgeLocalID, :) = EdgeAttr(EdgeAcc(rowIndex, 5), :);
            %EdgeAttr(GlobalEdgeID, 3:5) = [i, j, edgeLocalID]; 
            continue
        end
        
        % boundary conditions
        if (points(ID1, 3) == 0 && points(ID2, 3) == 0)
            EdgeAttr(edgeLocalID, 2) = 2;
            EdgeAttr(edgeLocalID, 3) = 0;
            NumDirichiletEdge = NumDirichiletEdge + 1;
            continue;
        end

        if (points(ID1, 3) == 1 && points(ID2, 3) == 1)
            EdgeAttr(edgeLocalID, 2) = 2;
            EdgeAttr(edgeLocalID, 3) = 1;
            NumDirichiletEdge = NumDirichiletEdge + 1;
            continue;
        end
        
        if (points(ID1, 1) == 0 && points(ID2, 1) == 0) || ...
           (points(ID1, 1) == 1 && points(ID2, 1) == 1) || ...
           (points(ID1, 2) == 0 && points(ID2, 2) == 0) || ...
           (points(ID1, 2) == 1 && points(ID2, 2) == 1)
            EdgeAttr(edgeLocalID, 2) = 3;
            EdgeAttr(edgeLocalID, 3) = 0.;
            NumNeumannEdge = NumNeumannEdge + 1;
            continue;
        end

        % common interface edge
        EdgeAttr(edgeLocalID, 2) = 0;
        EdgeAttr(edgeLocalID, 3) = 0.;
        NumInteriosEdge = NumInteriosEdge + 1;
    end
end
GlobalEdgeID = GlobalEdgeID - 1;
NumInteriosEdge = NumInteriosEdge - 1;
NumNeumannEdge = NumNeumannEdge - 1;
NumDirichiletEdge = NumDirichiletEdge - 1;

inx = find(EdgeAttr(:, 2) == 2);
inx = ceil(inx ./ 3);
figure(6)
view(3)
title('D releted triangles')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', Tri_Frac(inx, 1:3), 'FaceVertexCData', zeros(size(inx, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
pbaspect([1, 1, 1])

inx = find(EdgeAttr(:, 2) == 3);
inx = ceil(inx ./ 3);
figure(7)
view(3)
title('N releted triangles')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
patch('Vertices', points, 'Faces', Tri_Frac(inx, 1:3), 'FaceVertexCData', zeros(size(inx, 1), 1), 'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on
pbaspect([1, 1, 1])

%--------------------------------
%--------------------------------
%--------------------------------
% q3, p3, q2, p2, p1

Dim = globalTriID + NumTets + NumFracTri * 3 + NumFracTri + 1 + NumInteriosEdge;
