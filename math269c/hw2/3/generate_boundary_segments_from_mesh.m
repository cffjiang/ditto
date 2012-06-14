function boundary_segments = generate_boundary_segments_from_mesh(elements,nodes)

boundary_segments = [];

N_nodes = size(nodes,1);
N_elements = size(elements,1);

edge_table = sparse(N_nodes,N_nodes); % edge_table(i,j) = 1 means edge i-j exists

for t = 1:N_elements
    edge1 = [elements(t,1) elements(t,2)];
    edge2 = [elements(t,2) elements(t,3)];
    edge3 = [elements(t,3) elements(t,1)];
    
    edge_table(edge1(1),edge1(2)) = 1;
    edge_table(edge2(1),edge2(2)) = 1;
    edge_table(edge3(1),edge3(2)) = 1;
end

[is,js,vs] = find(edge_table);

for e = 1:nnz(edge_table) % loop over non-zero entries
    i = is(e);
    j = js(e);
    v = vs(e);
    if v~=1
        display('Fatal error!')
    end
    if edge_table(j,i) == 0
        boundary_segments = [boundary_segments; i,j];
    end
end

end