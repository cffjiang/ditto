function [nodal_normals boundary_segment_normals] = compute_normals(elements,nodes,boundary_segments)

N_nodes = size(nodes,1);
N_elements = size(elements,1);
N_boundary_segments = size(boundary_segments,1);

nodal_normals = zeros(N_nodes,2);
boundary_segment_normals = zeros(N_boundary_segments,2);

for i = 1:N_boundary_segments
    node_a_index = boundary_segments(i,1);
    node_b_index = boundary_segments(i,2);
    
    segment_vec = nodes(node_b_index,:)-nodes(node_a_index,:);
    segment_normal = left_handed_perpendicular_vec_normalized(segment_vec);
    
    boundary_segment_normals(i,:) = segment_normal;
    
    nodal_normals(node_a_index,:) = nodal_normals(node_a_index,:)+0.5*segment_normal;
    nodal_normals(node_b_index,:) = nodal_normals(node_b_index,:)+0.5*segment_normal;
end

for i = 1:N_nodes
    if nodal_normals(i,:)~=[0 0]
        nodal_normals(i,:) = nodal_normals(i,:)/norm(nodal_normals(i,:));
    end
end

end


function result = left_handed_perpendicular_vec_normalized(v)

    result = [v(2), -v(1)];
    result = result/norm(result);

end