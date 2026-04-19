function [mesh]=meshstruct2(msh)

        

        p=msh.POS';
        t=msh.TRIANGLES';
        mesh.p = p;
        mesh.t = t;
        mesh.np = length(unique(mesh.t(1:end-1,:)));
        mesh.ne = length(t); 
        mesh.nt = mesh.ne;
        mesh.N = 6; %Dof per node

end