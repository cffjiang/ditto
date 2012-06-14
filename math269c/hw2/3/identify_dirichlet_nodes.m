function [dirichlet_data dirichlet_node_list]= identify_dirichlet_nodes(nodes, xmin, xmax, ymin, ymax, dirichlet_value)
% dirichlet nodes are nodes cutted out from the input box
% dirichlet_data:
% the first column is bool, flag of whether a node is dirichlet
% the second column is the corresponding dirichlet value

N = size(nodes,1);

dirichlet_data = zeros(N,2);
dirichlet_node_list = [];

for i = 1:N
    x = nodes(i,1);
    y = nodes(i,2);
    
    if x<xmin || x>xmax || y<ymin || y>ymax
        dirichlet_data(i,1) = 1;
        dirichlet_data(i,2) = dirichlet_value;
        
        dirichlet_node_list = [dirichlet_node_list, i];
    end
end

end