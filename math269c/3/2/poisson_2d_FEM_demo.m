%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Poisson 2d FEM solver
%
% Solves -div(grad(u)) = f
%        n dot grad(u) = g (full Neumann boundary)
% over a triangle mesh
%
% Equivalent to solving Au=G+F, where Aij = int_Omega{Grad_Ni dot Grad_Nj}
%                                     Gi  = int_PartialOmega{g Ni}
%                                     Fi  = int_Omega{f Ni}
%
% TODO: let it support dirichlet boundary condition and mixed boundary condition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; close all; clc
fprintf('###  Poisson 2D FEM Solver with full Neumann boundary condition  ###\n')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% input the problem
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_number = 3;

if test_number == 1
    u_exact = @(x,y) x+2*y;                   % exact solution
    grad_u_exact = @(x,y) [1, 2];       % exact grad(u)
    f = @(x,y) 0;                      % exact forcing term
    
elseif test_number == 2
    u_exact = @(x,y) x^2+y^2;                   % exact solution
    grad_u_exact = @(x,y) [2*x, 2*y];       % exact grad(u)
    f = @(x,y) -4;                      % exact forcing term

elseif test_number == 3
    u_exact = @(x,y) x^3+y^3+x^2+y^2;                   % exact solution
    grad_u_exact = @(x,y) [3*x^2+2*x, 3*y^2+2*y];       % exact grad(u)
    f = @(x,y) -(6*x+2+6*y+2);                      % exact forcing term
end

fprintf('The exact solution is: \n')
display(u_exact)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% some pre-processing for the mesh
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tic;

elements = load('mesh_with_holes.dat');                                      % element data
N_elements = size(elements,1);
nodes = load('nodes.dat');                                                   % node data
N_nodes = size(nodes,1);
time_cost = toc; display(strcat('Time cost for loading mesh: ',num2str(time_cost))); tic;

boundary_segments = generate_boundary_segments_from_mesh(elements,nodes);    % boundary segment data
boundary_nodes = boundary_segments(:,1);                                     % boundary node data
time_cost = toc; display(strcat('Time cost for computing boundary segment mesh: ',num2str(time_cost))); tic;

[nodal_normals boundary_segment_normals] = compute_normals(elements,nodes,boundary_segments); % outward normals on nodes (zero for non-boundary nodes) and on boundary segs
time_cost = toc; display(strcat('Time cost for computing boundary segment normals: ',num2str(time_cost))); tic;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% find the proper input Neumann boundary condition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[g_nodal g_boundary] = compute_flux(elements,nodes,boundary_segments,nodal_normals,boundary_segment_normals,grad_u_exact);    % g is defined on all nodes, but we only care about those on boundary
time_cost = toc; display(strcat('Time cost for computing Neumann boundary condition: ',num2str(time_cost))); tic;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% build and solve FEM system Au=G+F
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[A G F] = build_system(elements,nodes,boundary_segments,f,g_boundary);
rhs = G+F;
time_cost = toc; display(strcat('Time cost for building linear system: ',num2str(time_cost))); tic;

u = minres(A,rhs,1e-10,1000);
time_cost = toc; display(strcat('Time cost for solving with minres: ',num2str(time_cost))); tic;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% plot the final solution u
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% compute the 3d plot coordinates
for i = 1:size(nodes,1)
    x(i) = nodes(i,1);
    y(i) = nodes(i,2);
    z(i) = u(i);
    z_exact(i) = u_exact(nodes(i,1),nodes(i,2));
end

% shift u and u_exact by their minimum: because if u satisfy the equation, then u+C also does
min_z = min(z);
min_z_exact = min(z_exact);
for i = 1:size(nodes,1)
    z(i) = z(i)-min_z;
    z_exact(i) = z_exact(i)-min_z_exact;
end

% plot the material space
figure
triplot(elements,x,y);
title('material space mesh')

% plot numerical solution
figure
subplot(1,2,1)
trisurf(elements,x,y,z,'EdgeColor','none');
title('numerical solution')

% plot exact solution
subplot(1,2,2)
trisurf(elements,x,y,z_exact,'EdgeColor','none');
title('exact solution')

time_cost = toc; display(strcat('Time cost for plotting: ',num2str(time_cost))); tic;

fprintf('Error: max_abs(u - u_exact) = %f\n', max(abs(z - z_exact)))
fprintf('END\n')