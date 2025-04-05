% The book is Numerical Analysis by Richard L. Burden and J. Douglas Faires
% This is Example 2 of Burden's book, page 722
% The exact solution is 
% u(x,t)= x*exp(y)
clear
clc
% Initializing.
nx = 7;            % Number of interior points in the x dimension.
ny = 6;            % Number of interior points in the y dimension.
Lx = 2;            % Length of the domain in the x direction.
Ly = 1;            % Length of the domain in the y direction.
hx = Lx / (nx-1);  % Grid spacing in x.
hy = Ly / (ny-1);  % Grid spacing in y.
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);

% Initialize A and b.
N = nx * ny;       % Total number of unknowns for non-square grid.
A = zeros(N, N);
b = zeros(N, 1);

% Source term f(x, y).
f = @(x, y) x*exp(y); 

% Boundary conditions.
bc_left = 0;       % u(x=0, y).
bc_right = 2;      % u(x=Lx, y).
bc_bottom = 1;     % u(x, y=0).
bc_top = exp(1);   % u(x, y=Ly).

% Index mapping function.
eta = @(i, j) j + (i - 1) * ny;  % Maps 2D indices to 1D index (vertical).

% Construct the system of equations.
for i = 1:nx
    for j = 1:ny
        idx = eta(i, j);  % Map (i, j) to 1D index.
        % Boundary conditions.
        if i == 1 || i == nx  || j == 1 || j == ny
            A(idx, idx) = 1;
            if i == 1
                b(idx) = bc_left;
            elseif i == nx
                b(idx) = bc_right *exp( y(j));
            end
            if j == 1
                b(idx) = bc_bottom*x(i);
            elseif j == ny
                b(idx) = bc_top * x(i);
            end
        else
            % Interior points.
            A(idx, eta(i-1, j)) = 1 / hx^2;
            A(idx, eta(i, j-1)) = 1 / hy^2;
            A(idx, idx) = -2 * (1 / hx^2 + 1 / hy^2);
            A(idx, eta(i, j+1)) = 1 / hy^2;
            A(idx, eta(i+1, j)) = 1 / hx^2;
            % Add source term to b.
            b(idx) = f(x(i), y(j));
        end
    end
end

% Solve the system.
w = A \ b;

% Reshape solution to 2D grid.
W = reshape(w, [ny, nx]);

% Display the solution at each node.
disp('Solution at each node:');
for i = 1:nx
    for j = 1:ny
        fprintf('node(%d) = (x(%d), y(%d)) = (%.2f, %.2f) = %.2f \n', ...
            eta(i, j), i, j, x(i), y(j), W(j, i));
    end
end

% Optional.
% Display the result
%disp('The solution at the grid points is:');
%disp(W);

% Exact Solution with error table.
sol=@(x,y) x*exp(y); 
exact_sol=zeros(N,1);
for i=1:nx
    for j=1:ny
        n=j+(i-1)*ny;
        exact_sol(n)=sol(x(i),y(j));
    end
end
err=abs(exact_sol-w);
table(w,exact_sol,err)

% Contour Plot.
Exact_sol=reshape(exact_sol,[ny,nx]);
figure()
contourf(W)
title('Numerical')
figure()
contourf(Exact_sol)
title('Exact')

% Surface Plot.
 figure()
 [X, Y] = meshgrid(x, y);
 surf(X, Y, W);
 xlabel('x');
 ylabel('y');
 zlabel('u(x,t)');
 title('Wave Equation Solution by Finite Difference Scheme','fontsize',14)
