
 % Solving: u_tt = u_xx, 0 < x < 1.
 % IC: u(x,0) = sin(pi*x), u_t(x,0) = cos(pi*x).
 % Neumann BC: u_x(0,t) = g(t), u_x(1,t) = h(t).
 clear
 clc
 % Initializing
 L = 1;                % Length of the domain.
 T = 1;                % Total time.
 nx = 51;              % Number of spatial points.
 nt = 101;             % Number of time points.
 hx = L / nx;          % Spatial step size.
 ht = T / nt;          % Time step size.
 lambda = (ht / hx)^2; % CFL number.
 
 % Stability check
 if lambda > 1
 error('Stability condition violated: lambda must be <= 1');
 end
% We also can use the Von Neumann condition where lambda must be less 
% than 0.5.
% Construct time and space vectors.

 x = linspace(0, L, nx);
 t = linspace(0, T, nt);
 % Initial and boundary conditions.
 
 f = @(x) sin(pi * x);  % u(x, 0).
 g = @(t) 1;            % Neumann BC: u_x(0,t) = g(t).
 h = @(t) 0;            % Neumann BC: u_x(1,t) = h(t).
 v = @(x) cos(pi * x);  % Initial velocity u_t(x, 0).
 
 % Initialize solution matrix.
 U = zeros(nx, nt);
 
 % 1- Apply initial condition u(x,0).
 U(:, 1) = f(x);
 
 % 2- Apply boundary conditions for the first time step.
 U(1, 1)  = U(2, 1)- hx * g(t(1));      % Neumann BC at x = 0.
 U(nx, 1) = U(nx-1, 1) + hx * h(t(1));  % Neumann BC at x = 1.
 
 % 3- Compute u(x, dt) using the given initial velocity.
 for i = 2:nx-1 %better approximate of u_t(x,0).
 U(i, 2) = (1- lambda) * f(x(i)) + ...
 0.5 * lambda * (f(x(i+1)) + f(x(i-1))) + ...
 ht * v(x(i));
 end
 
 % Enforce Neumann boundary conditions at second time step.
 U(1, 2) = U(2, 2) - hx * g(t(2));       % Neumann BC at x = 0.
 U(nx, 2) = U(nx-1, 2) + hx * h(t(2));   % Neumann BC at x = 1.
 
 % 4- Time-stepping for u(x, t) using finite-difference method.
 for j = 2:nt-1
 for i = 2:nx-1
 U(i, j+1) = 2 * (1- lambda) * U(i, j) + ...
 lambda * (U(i+1, j) + U(i-1, j))- ...
 U(i, j-1);
 end
 % Enforce Neumann boundary conditions at each time step.
 U(1, j+1) = U(2, j+1)- hx * g(t(j+1));        % Neumann BC at x = 0.
 U(end, j+1) = U(end-1, j+1) + hx * h(t(j+1)); % Neumann BC at x = 1.
 end
 
 % Surface Plot.
 [X, T] = meshgrid(x, t);
 surf(X, T, U');
 xlabel('x');
 ylabel('t');
 zlabel('u(x,t)');
 title('Wave Equation Solution with Neumann Boundary Conditions');
 