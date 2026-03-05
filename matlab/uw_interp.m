function []=uw_interp()
%UW_INTERP Interpolate grid using nearest neighbour (HPC Optimized)
%
%   Usage: 
%       uw_interp()
%       (Reads 'uw_grid.mat' and outputs 'uw_interp.mat')
%
%   ======================================================================
%   MODIFICATION HEADER (StaMPS-HPC)
%   ======================================================================
%   Author:        Mingjia Li
%   Date:          February 2026
%   Version:       1.0 
%   License:       GPL v3.0 (Inherited from StaMPS)
%
%   HPC Optimization:
%   1. Algorithm Replacement: Replaced external system call 'triangle' 
%      with MATLAB native 'delaunayTriangulation' to remove OS dependency 
%      and improve portability.
%   2. Search Efficiency: Replaced legacy 'dsearchn' with object-based 
%      'nearestNeighbor' method for faster and more robust lookups.
%
%   ======================================================================
%   ORIGINAL HEADER (StaMPS)
%   ======================================================================
%   Original Author: Andy Hooper, May 2007
%   ======================================================================

fprintf('Interpolating grid (HPC Version)...\n')

% 1. Load Data
if exist('uw_grid.mat', 'file')
    uw = load('uw_grid', 'n_ps', 'nzix');
else
    error('uw_grid.mat not found.');
end

[nrow, ncol] = size(uw.nzix);
[y, x] = find(uw.nzix);

% 2. Build Triangulation
% Construct Delaunay triangulation from PS pixel coordinates
DT = delaunayTriangulation(double(x), double(y)); 

% 3. Nearest Neighbor Search
% Create full grid query points
[X, Y] = meshgrid(1:ncol, 1:nrow);
query_points = [X(:), Y(:)];

% Perform search using optimized object method
Zvec = nearestNeighbor(DT, query_points);

% 4. Construct Grid Topology
% Define vertical (column) edges
grid_edges = [Zvec(1:end-nrow), Zvec(nrow+1:end)]; 

% Define horizontal (row) edges
Z_transposed_vec = reshape(reshape(Zvec, nrow, ncol)', nrow*ncol, 1);
grid_edges = [grid_edges; [Z_transposed_vec(1:end-ncol), Z_transposed_vec(ncol+1:end)]]; 

% Sort edges to canonical order (lowest node index first)
[sort_edges, I_sort] = sort(grid_edges, 2); 
edge_sign = I_sort(:,2) - I_sort(:,1);

% Remove duplicate edges and self-loops
[alledges, I, J] = unique(sort_edges, 'rows'); 
sameix = (alledges(:,1) == alledges(:,2));
alledges(sameix, :) = 0; % Mark self-loops

% Finalize unique edge list (dropping 0-index if present)
[edgs, I2, J2] = unique(alledges, 'rows');
n_edge = size(edgs, 1) - 1;
edgs = [[1:n_edge]', edgs(2:end, :)]; 

% Map edges back to grid indices (rowix/colix)
gridedgeix = (J2(J) - 1) .* edge_sign; 
colix = reshape(gridedgeix(1:nrow*(ncol-1)), nrow, ncol-1);
rowix = reshape(gridedgeix(nrow*(ncol-1)+1:end), ncol, nrow-1)';

fprintf('   Number of unique edges in grid: %d\n', n_edge);

% 5. Save Output
Z = reshape(Zvec, nrow, ncol);
save('uw_interp','edgs','n_edge','rowix','colix','Z');

end