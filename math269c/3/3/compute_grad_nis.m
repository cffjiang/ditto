function grad_nis = compute_grad_nis(node_pos1,node_pos2,node_pos3,element_area)

    grad_nis = zeros(3,2);
    twice_area = 2*element_area;
    
    grad_nis(1,1) = (node_pos2(2)-node_pos3(2)) / twice_area;
    grad_nis(1,2) = (node_pos3(1)-node_pos2(1)) / twice_area;
    grad_nis(2,1) = (node_pos3(2)-node_pos1(2)) / twice_area;
    grad_nis(2,2) = (node_pos1(1)-node_pos3(1)) / twice_area;
    grad_nis(3,1) = (node_pos1(2)-node_pos2(2)) / twice_area;
    grad_nis(3,2) = (node_pos2(1)-node_pos1(1)) / twice_area;
    
end