function [g_nodal g_boundary] = compute_flux(elements,nodes,boundary_segments,nodal_normals,boundary_segment_normals,grad_u_exact)

N_nodes = size(nodes,1);
N_elements = size(elements,1);
N_boundary_segments = size(boundary_segments,1);

g_nodal = zeros(N_nodes,1);
g_boundary = zeros(N_boundary_segments,1);

for b = 1:N_boundary_segments
    node_index = boundary_segments(b,1);
    node_position = nodes(node_index,:);
    x = node_position(1);
    y = node_position(2);
    g_nodal(node_index) = dot(grad_u_exact(x,y), nodal_normals(node_index,:));
    
    boundary_mid_point = 0.5 * (nodes(boundary_segments(b,1),:) + nodes(boundary_segments(b,2),:));
    x = boundary_mid_point(1);
    y = boundary_mid_point(2);
    g_boundary(b) = dot(grad_u_exact(x,y), boundary_segment_normals(b,:));
end

end