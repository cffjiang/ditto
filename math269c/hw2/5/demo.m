clear all; close all; clc

% forcing term
f = @(x) 24;

% boundary conditions
u0 = 1;
g = 14;

% non-uniform dx
N = 20;
dx_average = 1/(N-1);
dx = ones(N-1,1)*dx_average;

% initialize mesh
x = zeros(1,N);
x(1) = 0;
for i = 2:N
    x(i) = sum(dx(1:i-1));
end

% exact solution
u = zeros(1,N);
for i = 1:N
    u(i) = 12*x(i)^2  - 10*x(i) + 1;
end

% initial guess
v = zeros(1,N);
v(1) = u0;

% build rhs
rhs = zeros(N-1,1);
for i = 1:N-2
    rhs(i) = -f(x(i+1))*dx_average;
end
rhs(N-1) = -f(x(N))*dx_average/2+ g;
rhs(1) = rhs(1) + u0/dx_average; % contribution from dirichlet

% build matrix
A = sparse(N-1,N-1);
A(1,1) = A(1,1) + 1/dx_average; % contribution from dirichlet
for e = 2:N-1
    dx_current = dx_average;
    A(e-1,e-1) = A(e-1,e-1) + 1/dx_current;
    A(e-1,e) = A(e-1,e) - 1/dx_current;
    A(e,e-1) = A(e,e-1) - 1/dx_current;
    A(e,e) = A(e,e) + 1/dx_current;
end

% solve for solution
unknown = A\rhs;

% build v
v(2:N) = unknown;

% plot
figure()
plot(x,u,'.-',x,v,'o-');

max(abs(u-v))




