function [mesh]=meshstruct(mesh,p,e,t,g);

        mesh.p = p;
        mesh.t = t;
        mesh.np = length(unique(mesh.t(1:end-1,:)));
        mesh.p = [mesh.p;zeros(1,mesh.np)];
        mesh.ne = length(t);    
        mesh.np = length(unique(t(1:end-1,:)));
        mesh.nt = mesh.ne;
        mesh.N = 6; %Dof per node
        mesh.e=e;
        mesh.g=g;

%         p=msh.POS';
%         t=msh.TRIANGLES';
%         mesh.p = p;
%         mesh.t = t;
%         mesh.np = length(unique(mesh.t(1:end-1,:)));
%         mesh.ne = length(t); 
%         mesh.nt = mesh.ne;
%         mesh.N = 6; %Dof per node

end