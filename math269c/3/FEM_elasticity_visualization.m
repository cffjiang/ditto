clear
mesh=load('mesh_with_holes.dat');
vertices=load('nodes.dat');

s=size(vertices);
num_verts=s(1);
mesh_size=size(mesh);

u=zeros(2*size(vertices,1),1);
for i=1:num_verts
    u(2*i) = sin(vertices(i,1));
end

stress=zeros(mesh_size(1),4);
for t=1:mesh_size
    % if vertices(mesh(t,1),1) < 0.5
        % stress(t,1)=2;
    % else
        % stress(t,1)=1;
    % end
end




node_stress=zeros(num_verts,1);

for t=1:mesh_size(1)
    s_norm=sqrt(stress(t,1)*stress(t,1)+stress(t,2)*stress(t,2)+stress(t,3)*stress(t,3)+stress(t,4)*stress(t,4));
    for i=1:3
        node_stress(mesh(t,i))=node_stress(mesh(t,i))+s_norm/3;
    end
end

for i=1:num_verts
    x(i)=vertices(i,1);
    y(i)=vertices(i,2);
    
    x_deformed(i)=x(i)+u(2*i-1);
    y_deformed(i)=y(i)+u(2*i);
end

figure

subplot(1,2,1)
triplot(mesh,x,y);

subplot(1,2,2)

trisurf(mesh,x_deformed,y_deformed,node_stress,'EdgeColor','none');
view(2)