%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Linear elasticity equilibrium 2d FEM solver
%
% Solves div(sigma) + f = 0
%        u = 0 at dirichelt boundary
% over a triangle mesh
%
% TODO: let it support arbitrary dirichlet value (currently only support 0)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; close all; clc
fprintf('###  Linear elasticity equilibrium 2D FEM Solver with zero Dirichlet boundary condition  ###\n')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% input the problem
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

youngs_modulus = 100;
poisson_ratio = 0.3;
density = 1;

dirichlet_box = [1e-8, 999, -999, 999]; % this is a [xmin xmax ymin ymax] box that cuts out dirichlet nodes. 
dirichlet_value = 0; % currently this program only supports 0.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% material and world parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lambda = youngs_modulus*poisson_ratio / ((1 + poisson_ratio)*(1 - 2*poisson_ratio));
mu = youngs_modulus / ( 2 * (1 + poisson_ratio) );
gravity = density * [0, -9.8]';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% some pre-processing for the mesh
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tic;

% load mesh from file
elements = load('mesh_with_holes.dat');                                      % element data
N_elements = size(elements,1);
nodes = load('nodes.dat');                                                   % node data
N_nodes = size(nodes,1);
time_cost = toc; display(strcat('Time cost for loading mesh: ',num2str(time_cost))); tic;

% compute boundary segment mesh
boundary_segments = generate_boundary_segments_from_mesh(elements,nodes);    % boundary segment data
boundary_nodes = boundary_segments(:,1);                                     % boundary node data
time_cost = toc; display(strcat('Time cost for computing boundary segment mesh: ',num2str(time_cost))); tic;

% identify dirichlet boundary
[dirichlet_data dirichlet_node_list] = identify_dirichlet_nodes(nodes, dirichlet_box(1), dirichlet_box(2), dirichlet_box(3), dirichlet_box(4), dirichlet_value);
time_cost = toc; display(strcat('Time cost for identifying dirichlet nodes: ',num2str(time_cost))); tic;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% build and solve FEM system
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[K rhs] = build_system(elements,nodes,dirichlet_data,dirichlet_node_list,lambda,mu,gravity);
time_cost = toc; display(strcat('Time cost for building linear system: ',num2str(time_cost))); tic;

u = minres(K,rhs,1e-10,1000);
time_cost = toc; display(strcat('Time cost for solving linear system: ',num2str(time_cost))); tic;

stress = evaluate_stress(elements,nodes,u,lambda,mu);
time_cost = toc; display(strcat('Time cost for evaluating elementwise Cauchy stress: ',num2str(time_cost))); tic;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% plot the final solution u
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% compute the node stress
node_stress=zeros(N_nodes,1);
for t=1:N_elements
    s_norm=sqrt(stress(t,1)*stress(t,1)+stress(t,2)*stress(t,2)+stress(t,3)*stress(t,3)+stress(t,4)*stress(t,4));
    for i=1:3
        node_stress(elements(t,i))=node_stress(elements(t,i))+s_norm/3;
    end
end

% compute world space coordinates
for i=1:N_nodes
    x(i)=nodes(i,1);
    y(i)=nodes(i,2);
    x_deformed(i)=x(i)+u(2*i-1);
    y_deformed(i)=y(i)+u(2*i);
end

% plot the material space
figure
triplot(elements,x,y);
title('material space mesh')

% plot the world space with stress color
figure
trisurf(elements,x_deformed,y_deformed,node_stress,'EdgeColor','none');
view(2)

time_cost = toc; display(strcat('Time cost for plotting: ',num2str(time_cost))); tic;
fprintf('END\n')