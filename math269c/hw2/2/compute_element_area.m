function area = compute_element_area(elements,nodes,e)

    index_a = elements(e,1);
    index_b = elements(e,2);
    index_c = elements(e,3);

    A = nodes(index_a,:);
    B = nodes(index_b,:);
    C = nodes(index_c,:);
    
    vec1 = B-A;
    vec2 = C-A;
    
    area = 0.5 * abs(vec1(1)*vec2(2)-vec1(2)*vec2(1));

end