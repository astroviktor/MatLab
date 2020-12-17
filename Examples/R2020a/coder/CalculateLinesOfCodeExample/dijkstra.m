% DIJKSTRA Find length of shortest path between nodes in a graph
%
% D = dijkstra(A, p) 
% Takes a graph represented by its adjacency matrix 'A' along with a node 
% 'p' as input and returns a vector 'D' containing the length of the 
% shortest path from 'p' to all other nodes in the graph. 

% Copyright 2018 The MathWorks, Inc.
function D = dijkstra(A, p) %#codegen

    narginchk(2,2);
    
    [m, n] = size(A);

    % Assertions to make sure inputs are valid
    assert(m == n, "Input adjacency matrix for graph must be a square matrix");
    assert(rem(p, 1) == 0 && p <= m && p > 0, "Input src must be a node in the graph");
    
    % Initialization
    max = realmax;
    D = repmat(max, 1, m);
    D(p) = 0;
    visited = false(1, m);
    
    for i = 1:m
        % Select next node to visit
        min = max;
        u = -1;
        for v = 1:n
            if ~visited(v) && D(v) <= min
                min = D(v);
                u = v;
            end
        end
        
        % Mark selected node as visited
        visited(u) = true;
        
        %{ 
          Update distances of nodes adjacent to selected node that are yet
          to be visited
        %}
        for v = 1:n
            if(~visited(v) && A(u, v) ~= 0 && D(u) ~= max)
                distVal = D(u) + A(u, v);
                if distVal < D(v)
                    D(v) = distVal;
                end
            end
        end
    end
end