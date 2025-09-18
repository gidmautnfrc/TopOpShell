function [mesh, params, psi0, bc, signatures, pdecoef, matprop] = casquete(mesh,psi,params)
    %%cargar malla
   
mallas = {'casquete_mallag.m','casquete_mallam.m','casquete_mallaf.m','casquete_mallaf02.m','casquete_mallaf03.m'};

   if isempty(mesh)
	
	    index = 4;
        run(mallas{index});
		
        t = msh.TRIANGLES';
		
        psi0=-ones(length(unique(t(1:end-1,:))),1);

    else
		index = mesh.index;
		index = index + 1;
		
        run(mallas{index});

		[psi]=InterpolationCS(psi,mesh,msh);
		
        psi0 = psi;
    end

    %% Mesh struct
        
        [mesh]=meshstruct2(msh);
        
	    mesh.index = index;

    %rotaciones de elementos

        [mesh]=rotation(mesh);
    
    %% Material properties
  
    E = 70.0E+9; 
    gamma = 1.0e-4;
    nu = 0.3;
    h = 0.003; 
    
    matprop.E0 = E;
    matprop.gamma = gamma;
    matprop.nu0 = nu;
    matprop.h0 = h;
    matprop.E1 = gamma*E;
    
    
    %DKT (Discrete Kirchhoff Triangle)
    
    la = nu*E/((1+nu)*(1-2*nu)); 
    mu = E/(2*(1+nu)); % plane strain
    la = 2*mu*la/(la+2*mu); % plane stress
    
    matprop.la_dkt = la; 
    matprop.mu_dkt = mu;
    
    %%ANDES 
    matprop.la_andes = ((nu.*E.^2)./((1+nu).^2).*(1-2.*nu))./(((nu.*E)./(1+nu).*(1-2.*nu))+(E./(1+nu)));
    matprop.mu_andes = E./(2.*(1+nu));
    
    %% Coeficiente c pde
    
    [pdecoef] = cprop(matprop);    
    
    %%Firma andes   

    [signatures]=andes_signature('OPT',nu);
   
    %% parameters
        % line-search procedure
        params.kmin = 1.0E-3;
        % stop criterion
        params.stop = 1.0*pi/180;
        % volume penalization
        params.penalty = 1250; 
    %% boundary conditions
    
    bc.pNeu = []; bc.pDir = []; 
    x = mesh.p(1,:); y = mesh.p(2,:); z = mesh.p(3,:);
    eps = 0.01*sqrt(min(mesh.A));

    % XY Plane

         cz = 0; 
         node_z = find(z < cz + eps); node_z = node_z'; %Defino un vector columna con los indices de los nodos a restringir.
		 npb_z = size(node_z,1); node_z = [node_z;node_z;node_z];
         dof_z = [3*ones(npb_z,1);4*ones(npb_z,1);5*ones(npb_z,1)];
		 val_z = [zeros(npb_z,1);zeros(npb_z,1);zeros(npb_z,1)];
		 bc_z = [node_z, dof_z, val_z];
	
    % XZ Plane

		 cy = 0; 
         node_y = find(y < cy + eps); node_y = node_y';
		 npb_y = size(node_y,1); node_y = [node_y;node_y;node_y];
         dof_y = [2*ones(npb_y,1);4*ones(npb_y,1);6*ones(npb_y,1)];
		 val_y = [zeros(npb_y,1);zeros(npb_y,1);zeros(npb_y,1)];
		 bc_y = [node_y, dof_y, val_y];

    % ZY Plane	 

		 cx = 0; 
         node_x = find(x < cx + eps); node_x = node_x';
		 npb_x = size(node_x,1); node_x = [node_x;node_x;node_x];
         dof_x = [1*ones(npb_x,1);5*ones(npb_x,1);6*ones(npb_x,1)];
		 val_x = [zeros(npb_x,1);zeros(npb_x,1);zeros(npb_x,1)];
		 bc_x = [node_x, dof_x, val_x];
		 
		 bc.pDir = [ bc_z;
		             bc_y;
					 bc_x ];
        
%% cond de neumman

% Carga puntual aplicada en el centro del alma

   cx = 0.3535534; cy = 0.3535534; cz=0.7500000;
	
    [~,node1] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);

    bc.pNeu = [node1, 3, -150]; 
end