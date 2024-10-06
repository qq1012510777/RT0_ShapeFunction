function V = tetrahedronVolume(A, B, C, D)
    % Function to calculate the volume of a tetrahedron given four vertices
    % A, B, C, D are 1x3 vectors representing the coordinates of the vertices
    % Example: A = [x1, y1, z1], B = [x2, y2, z2], etc.
    
    % Calculate vectors relative to vertex A
    AB = B - A;
    AC = C - A;
    AD = D - A;
    
    % Create a matrix with AB, AC, AD as rows
    M = [AB; AC; AD];
    
    % Calculate the determinant of the matrix
    detM = det(M);
    
    % Calculate the volume
    V = abs(detM) / 6;
end