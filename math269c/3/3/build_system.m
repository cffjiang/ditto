function [K rhs] = build_system(elements,nodes,dirichlet_data,dirichlet_node_list,lambda,mu,gravity)

N_elements = size(elements,1);
N_nodes = size(nodes,1);

K = sparse(2*N_nodes,2*N_nodes);
rhs = zeros(2*N_nodes,1);

% build matrix
for e = 1:N_elements
    element = elements(e,:);
    X1 = nodes(element(1),:);
    X2 = nodes(element(2),:);
    X3 = nodes(element(3),:);
    
    area = compute_element_area(elements,nodes,e);
    
    Dm = [(X2-X1)' (X3-X1)'];
    Dm_inv = inv(Dm);
    a = Dm_inv(1,1);
    b = Dm_inv(1,2);
    c = Dm_inv(2,1);
    d = Dm_inv(2,2);
    
    M = [-(a+c)   0   a    0   c    0;
         0     -(a+c) 0    a   0    c;
         -(b+d)   0   b    0   d    0;
         0     -(b+d) 0    b   0    d];

    M_hat = [1   0    0    0;
             0  0.5  0.5   0;
             0  0.5  0.5   0;
             0   0    0    1] * M;
    Ke = M_hat' * [2*mu+lambda  0    0      lambda ;
                       0       2*mu  0          0  ;
                       0        0   2*mu        0  ;
                     lambda     0    0     2*mu+lambda] * M_hat;
    
    for ie = 1:3
        i = elements(e,ie);
        for je = 1:3
            j = elements(e,je);
            for a = 1:2
                for b = 1:2
                    K(2*(i-1)+a, 2*(j-1)+b) = K(2*(i-1)+a, 2*(j-1)+b) + area*Ke(2*(ie-1)+a, 2*(je-1)+b);
                end
            end
        end
    end
end
    
% build rhs
for e = 1:N_elements
    
    area = compute_element_area(elements,nodes,e);
    
    Fe = (1/3)*area*[gravity;gravity;gravity];
    
    for ie = 1:3
        i = elements(e,ie);
        rhs(2*(i-1)+2) = rhs(2*(i-1)+2) + Fe(2*(ie-1)+2);
    end
end

% impose zero dirichlet boundary
for it = 1:size(dirichlet_node_list,2)
    i = dirichlet_node_list(it);
    K(2*(i-1)+1,:) = 0;
    K(:,2*(i-1)+1) = 0;
    K(2*(i-1)+1,2*(i-1)+1) = 1;
    rhs(2*(i-1)+1) = 0;
    
    K(2*(i-1)+2,:) = 0;
    K(:,2*(i-1)+2) = 0;
    K(2*(i-1)+2, 2*(i-1)+2) = 1;
    rhs(2*(i-1)+2) = 0;
end


end % end of function: build_system