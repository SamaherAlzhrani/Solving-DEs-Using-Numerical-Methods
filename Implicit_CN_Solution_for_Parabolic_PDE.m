
% The book is Numerical Analysis by Richard L. Burden and J. Douglas
% Faires.
% Crank Niclson Finite Difference Method.
% This code solve parabolic PDE of the form ut- alpha uxx=2.
% This is Exercise 16 , section 12.2, page738.
% This scheme is unconditionally stable.
clear
clc
% 1- Initilizing.
L = 1;              % Length of the spatial domain.
T = 1;              % Total time.
C = 1;              % Thermal Diffusivity constant.
nx = 11;            % Number of spatial points.
nt = 101;           % Number of time steps.
hx = L / (nx -1);   % Spatial step size.
ht = T /(nt-1);     % Time step size.
lambda = C * ht / hx^2;

% Construct time and space vectors.
x= 0:hx:L;          % Space vector.
t= 0:ht:T;          % Time vector.

% The source of the PDE.
% 2- Initial and boundary conditions the grid of points.
g=@(x,t) 2;
gvec=zeros(nx-2,1);
for i=1:nx-2
gvec(i)=g(x(i+1),t(2));% Construct the source function vector at time and space.
end
u=zeros(nx,nt);

% Boundary conditions.
u(1,:)=0; u(nx,:)=0; 

% Initial condition.
f=@(x) sin(pi*x) + x.*(1-x);
u(:,1)=f(x);

% 3- The implicit method.
A=zeros(nx-2,nx-2); % The linear system we solve.
b_stage=u(2:nx-1,1)+ gvec * ht;

for k=2:nt % How many system to solve.
    % 1- fill matrix A.
    A(1,[1,2])=[1+lambda,-lambda/2]; % First row.
       j=1; % To fill the three desired columns of matrix A.
      for i=2:nx-3
       A(i,[j,j+1,j+2])=[-lambda/2,1+lambda,-lambda/2]; % inside the matrix.
          j=j+1;
      end
    A(nx-2,[nx-3,nx-2])=[-lambda/2,1+lambda]; %Last row.
   
    % 2- Fill matrix B, then find b by B*bb.
    
    B(1,[1,2])=[1-lambda,lambda/2]; % First row.
       j=1; % To fill the three desired columns of matrix B.
      for i=2:nx-3
       B(i,[j,j+1,j+2])=[lambda/2,1-lambda,lambda/2]; % inside the matrix.
          j=j+1;
      end
    % 3- Solve the system.
    B(nx-2,[nx-3,nx-2])=[lambda/2,1-lambda]; %Last row.
    b=B*b_stage;
    z=A\b;
    u(2:nx-1,k)=z;
    % 4- Update the source vector of the linear system.
    j=3; % update time index.
    for i=1:nx-2
    gvec(i)=g(x(i+1),t(j));% Construct the source function vector at time and space.
    end
    b_stage=z + gvec * ht;
    j=j+1;
end

% Exact Solution
usol=zeros(nx,nt);
sol=@(x,t) exp(-pi^2*t)*sin(pi*x) + x.*(1-x);
% dispaly book table at t=0.5
disp('This is the result of this equation for t=0.5, ')
approximate=u(:,(0.5/0.01+1));
exact=sol(x,0.5)';
err=abs(approximate-exact);
x=x';
table(x,approximate,exact,err)

% Plot 
 [X, T] = meshgrid(x, t);
 surf(X, T, u');
 xlabel('x');
 ylabel('t');
 zlabel('u(x,t)');
 title('Heat Equation Solution by the Crank-Niclson Difference Scheme.','fontsize',13)
