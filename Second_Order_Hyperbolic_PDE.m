% The book is Numerical Analysis by Richard L. Burden and J. Douglas Faires
% Finite Difference Method for 1D Hyperbolic PDE.
% This code solve Hyberpolic PDE of the form u_tt-(alpha)^2 u_xx=0.
% This is Example1, section 12.3, page743.
% The exact solution is u(x,t)=sin(pi*x)*cos(2*pi*t).
clear
clc
% Initializing.
L = 1;                 % Length of the spatial domain.
T = 1;                 % Total time.
alpha_square = 4;      % Wave speed.
nx = 11;               % Number of spatial points.
nt = 21;               % Number of time points.
hx = L / (nx -1);      % spatial step size.
ht = T /(nt-1);        % Time step size.
lambda =( alpha_square * ht^2) / hx^2;

% Construct time and space vectors.
x= 0:hx:L; 
t= 0:ht:T;

% Von Neumann Stability check.
 if ht > hx / sqrt(alpha_square) * sqrt(2)
 error('Stability condition violated: ht <= hx/sqrt(alpha_square)*sqrt(2)');
 end

% The source of the PDE.
g=@(x,t) 0;
u=zeros(nx,nt);
% Initial condition.

f=@(x) sin(pi*x);
u(:,1)=f(x);
% Neumann condition in this case, u_t=0.
u(:,2)=f(x);
% Boundary conditions.
u(1,:)=0; u(nx,:)=0; 

% Explicit Scheme.
for j=2:nt-1
    for i=2:nx-1
   u(i,j+1)=lambda*u(i-1,j)+(2-2*lambda)*u(i,j)+lambda*u(i+1,j)-u(i,j-1);     
    end
end

% Exact Solution.
usol=zeros(nx,nt);
sol=@(x,t) sin(pi*x)*cos(2*pi*t);
for i=1:nx
    for j=1:nt
        usol(i,j)=sin(pi*x(i))*cos(2*pi*t(j));
    end
end

% For displaying table.
 disp(' i |  j |   x   |   t    |  uapp  |  exact  |  error')
 disp('------------------------------------------------------------')
for i=1:nx
 for j=1:nt
 usol(i,j) = sol(x(i),t(j));
 err=abs(usol(i,j)-u(i,j));
 fprintf('%2.0f | %2.0f | %4.3f | %4.3f  |%7.5f | %7.5f | %11.10f \n' , i, j, x(i), t(j), u(i,j), usol(i,j), err)
 end
end

% Contour Plot.
figure()
contourf(u)
colormap(jet(256))
title('Numerical')
figure()
contourf(usol)
colormap(jet(256))
title('Exact')

% Surface Plot.
 figure()
 [X, T] = meshgrid(x, t);
 surf(X, T, u');
 xlabel('x');
 ylabel('t');
 zlabel('u(x,t)');
 title('1D Hyperbolic wave Equation Solution by the Finite Difference Scheme. ','fontsize',10)
