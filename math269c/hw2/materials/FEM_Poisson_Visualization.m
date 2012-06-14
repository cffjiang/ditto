
clear

mesh=load('mesh_with_holes.dat');
vertices=load('nodes.dat');
s=size(vertices);
num_verts=s(1);

%%%%%%%%%%%%%%%%%%%%%%%%%
u=ones(size(vertices,1),1);
for i=1:num_verts
    y(i)=vertices(i,2);
    u(i) = y(i)^2-y(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:num_verts
    x(i)=vertices(i,1);
    y(i)=vertices(i,2);
    z(i)=u(i);
end

figure

subplot(1,2,1)
triplot(mesh,x,y);

subplot(1,2,2)

trisurf(mesh,x,y,z,'EdgeColor','none');

% plot3(xs,ys,zs);