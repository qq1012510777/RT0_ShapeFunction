function points = generatePointsInTetrahedron(A, B, C, D, n)
    % Function to generate 'n' uniformly distributed points in a tetrahedron
    % A, B, C, D are 1x3 vectors representing the coordinates of the vertices.
    % 'n' is the number of points to generate.
    
    % Initialize an array to store the generated points
    points = zeros(n, 3);
    
    for i = 1:n
        % Generate 4 random numbers and sort them
        r = sort(rand(1, 3));
        
        % Compute the barycentric coordinates
        lambda1 = r(1);
        lambda2 = r(2) - r(1);
        lambda3 = r(3) - r(2);
        lambda4 = 1 - r(3);
        
        % Calculate the Cartesian coordinates of the point
        points(i, :) = lambda1 * A + lambda2 * B + lambda3 * C + lambda4 * D;
    end
end
