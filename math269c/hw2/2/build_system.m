%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% build_system
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [A G F] = build_system(elements,nodes,boundary_segments,f,g_boundary)

A = build_A(elements,nodes);
G = build_G(elements,nodes,boundary_segments,g_boundary);
F = build_F(elements,nodes,f);

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% build_A
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function A = build_A(elements,nodes)

N_elements = size(elements,1);
N_nodes = size(nodes,1);

A = sparse(N_nodes,N_nodes);

for e = 1:N_elements
    first_node = nodes(elements(e,1),:);
    second_node = nodes(elements(e,2),:);
    third_node = nodes(elements(e,3),:);
    
    element_area = compute_element_area(elements,nodes,e);
    
    grad_nis = compute_grad_nis(first_node,second_node,third_node,element_area);
    
    for a = 1:3
        for b = 1:3
            A(elements(e,a),elements(e,b)) = A(elements(e,a),elements(e,b)) + element_area*dot(grad_nis(a,:),grad_nis(b,:));
        end
    end
end

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% build_G
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function G = build_G(elements,nodes,boundary_segments,g_boundary)

N_elements = size(elements,1);
N_nodes = size(nodes,1);
N_boundary_segments = size(boundary_segments,1);

G = zeros(N_nodes,1);

for b = 1:N_boundary_segments
    first_node_index = boundary_segments(b,1);
    second_node_index = boundary_segments(b,2);
    first_node = nodes(first_node_index,:);
    second_node = nodes(second_node_index,:);
    
    segment_length = norm(second_node-first_node);
    
    G(first_node_index) = G(first_node_index) + 0.5*segment_length*g_boundary(b);
    G(second_node_index) = G(second_node_index) + 0.5*segment_length*g_boundary(b);
end

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% build_F
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function F = build_F(elements,nodes,f)

N_elements = size(elements,1);
N_nodes = size(nodes,1);

F = zeros(N_nodes,1);

for e = 1:N_elements
    first_node = nodes(elements(e,1),:);
    second_node = nodes(elements(e,2),:);
    third_node = nodes(elements(e,3),:);

    triangle_center_position = (first_node+second_node+third_node)/3;
    
    element_area = compute_element_area(elements,nodes,e);
    
    for a = 1:3
        F(elements(e,a)) = F(elements(e,a)) + (1/3)*element_area*f(triangle_center_position(1),triangle_center_position(2));
    end
end

end








