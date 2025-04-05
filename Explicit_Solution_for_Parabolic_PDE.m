clear 
clc
format short
% Example 1, page 727, of the book Numerical Analysis by Richard L. Burden
% and J. Douglas Faires
% ut-uxx=0, c=1, 0<x<1, t>=0
% hx=0.1, ht=0.0005
% The exact solution, u(x,t)=exp(-pi^2*t)*sin(pi*x)
disp('Forward-Difference method for solving Parabolic PDEs:')
disp('============================================================================')

% 1- Initilizing.

nx=11;               % The number of spatial points.
a=0; b=1;            % Spatila domain.
hx=(b-a)/(nx-1);     % The spatial step size.
x=a:hx:b;            % Spatial vector.
nt=2001;             % The number of spatial points.
t0=0; tend=1;        % Time domain.
ht=(tend-t0)/(nt-1); % The time step size.
t=t0:ht:tend;        % Time vector.
c=1;                 % Thermal diffusivity constant.
r=(c^2)*(ht/(hx)^2); % If r<0.5, then the explicit method is stable.

% Stability Check.
if r >= 0.5
    error('Stability condition violated: r must be < 0.5')
end

% 2- Initial and boundary conditions.

f=@(x)sin(pi*x);% initial condition f(x).
g1= @(t)0;      % boundary conditions g1(t) and g2(t).
g2= @(t)0;

% 3- Inilization of the grid of points
u= zeros(nx,nt);
u(:,1)= f(x);
u(1,:)=0; u(nx:1)=0;
% 4- The explicit method
 for j= 1:nt-1       % Time Loop 
     for i= 2:nx-1   % Spatial Loop
         u(i,j+1) = r*(u(i-1,j))+(1-2*r)*u(i,j) + r*u(i+1,j);
     end
 end

exact=zeros(nx,nt);
exact(:,1)= f(x);
exact(1,:)=0;exact(nx:1)=0;
for i=2:nx-1
    for j=2:nt
        exact(i,j)=(exp(-pi^2*t(j))*sin(pi*x(i)));
    end
end

% Plot the solution at t = 0.5.
[X, T] = meshgrid(0:hx:1, 0:ht:1);
contourf(X, T, u');
xlabel('x'); ylabel('t'); zlabel('u(x, t)');
title('Solution of 1D Unsteady Heat Equation');
plot(x',u(:,501))

% Displaying the result of the exact and approximated solution at t = 0.5.
%The errors
err1=abs(exact-u);
Approximate = u(:,1001);
Exact = exact(:,1001);
Error = err1(:,1001);
%table(x',Approximate,Exact,Error)
table=table(x',Approximate, Exact, Error, 'VariableNames', {'x', 'u(x,0.5)', 'Exact at t=0.5', 'Error'});
disp(table)