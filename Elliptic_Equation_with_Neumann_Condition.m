
clear
clc
% This code to solve Laplace equation uxx+uyy=0, Check the 
% related file.

% Parameters
nx = 3;             % Number of  points in x-direction.
ny = 5;             % Number of  points in y-direction.
Lx = 1;             % Length in x-direction.
Ly = 1;             % Length in y-direction.
hx = Lx / (nx - 1); % Step size in x.
hy = Ly / (ny - 1); % Step size in y.

% Coefficients.
cx = 1 / hx^2;           % Coefficient for x-direction.
cy = 1 / hy^2;           % Coefficient for y-direction.
c_diag = -2 * (cx + cy); % Diagonal coefficient.

% Initialize matrix A and vector b.
N = nx * ny;      % Total number of grid points.
A = sparse(N, N); % Sparse matrix for efficiency (similar to zeros matrix).
b = zeros(N, 1);  % Right-hand side vector.

% Construct the matrix A and vector b.
for k = 1:N
    % Simplified grid indices
    j = mod(k-1, ny) + 1; % Row index (y-direction).
    i = ceil(k / ny);     % Column index (x-direction).

    % Boundary conditions.
    if i == 1      % Left boundary (Dirichlet: u = 0)
        A(k, k) = 1;
        b(k) = 0;
    elseif i == nx % Right boundary (Dirichlet: u = 1)
        A(k, k) = 1;
        b(k) = 1;  % here condition
    elseif j == 1  % Bottom boundary (Neumann: du/dy = 0)
        A(k, k) = 1;
        A(k, k+1) = -1;
        b(k) = 0;
  %b(k) = -10*hy;  % in case Nuemann condition not zero and for example equal to 10.
    elseif j == ny % Top boundary (Dirichlet: u = 0)
        A(k, k) = 1;
        b(k) = 0;
    else % Interior points
        A(k, k) = c_diag; % Diagonal.
        A(k, k+1) = cy;   % Neighbor above.
        A(k, k-1) = cy;   % Neighbor below.
        A(k, k+ny) = cx;  % Neighbor to the right.
        A(k, k-ny) = cx;  % Neighbor to the left.
    end
end

% Solve the linear system.
u = A \ b;

% Reshape the solution vector into a 2D grid.
u_grid = reshape(u, ny, nx);

% Display the result (Optional).
%disp('Solution grid (u):');
%disp(u_grid');

% Plot the solution.
[X, Y] = meshgrid(0:hx:Lx, 0:hy:Ly);
surf(X, Y, u_grid);
xlabel('x'); ylabel('y'); zlabel('u(x, y)');
title('Solution to 2D Laplace Equation');
