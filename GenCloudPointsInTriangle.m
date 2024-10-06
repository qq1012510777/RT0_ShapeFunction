function [Points] = GenCloudPointsInTriangle(A, B, C, numPoints)
    % Generate random barycentric coordinates
    r1 = rand(numPoints, 1);
    r2 = rand(numPoints, 1);
    
    % Ensure the points are inside the triangle by taking square root
    r1_sqrt = sqrt(r1);
    u = 1 - r1_sqrt;
    v = r1_sqrt .* (1 - r2);
    w = r1_sqrt .* r2;
    
    % Compute the coordinates of the points in the triangle
    Points = u * A + v * B + w * C;
end