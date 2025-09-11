function [mesh]=meshstruct3(mesh)

        

        
        mesh.np = length(unique(mesh.t(1:end-1,:)));
        mesh.ne = length(mesh.t); 
        mesh.nt = mesh.ne;
        mesh.N = 6; %Dof per node

end