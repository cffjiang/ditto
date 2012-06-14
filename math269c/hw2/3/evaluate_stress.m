function stress = evaluate_stress(elements,nodes,u,lambda,mu)

N_elements = size(elements,1);
N_nodes = size(nodes,1);

stress = zeros(N_elements,4);

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

    u1 = [u(2*element(1)-1) u(2*element(1))];
    u2 = [u(2*element(2)-1) u(2*element(2))];
    u3 = [u(2*element(3)-1) u(2*element(3))];

    Du = [(u2-u1)' (u3-u1)'];
    
    du_dx = Du*Dm_inv;
    strain = 0.5*(du_dx+du_dx');
    sigma = 2*mu*strain + lambda*trace(strain);
    
    stress(e,:) = [sigma(1,1), sigma(1,2), sigma(2,1), sigma(2,2)];

end