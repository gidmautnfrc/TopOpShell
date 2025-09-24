function [mesh, params, psi0, bc, signatures, pdecoef, matprop] = paraboloid(mesh,psi,params)
    %%cargar malla
   
mallas = {'paraboloid_mallag.m','paraboloid_mallam.m','paraboloid_mallaf.m','paraboloid_mallaf02.m','paraboloid_mallaf03.m'};

   if isempty(mesh)
	
	    index = 1;
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
%         params.penalty = 200; 
        params.penalty = 6; 

    %% boundary conditions

    bc.pNeu = []; bc.pDir = []; 
    x = mesh.p(1,:); y = mesh.p(2,:); z = mesh.p(3,:);
    eps = 0.01*sqrt(min(mesh.A));

%    4 Apoyos fijos

%    cx = -0.5; cy = 0.5; cz=0;
% 
%     [~,node1] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
% 	node1 = [node1;node1;node1];
%     dof1 = [1;2;3];
% 	val1 = [0;0;0];	
% 	bc_1 = [node1, dof1, val1];
% 	
%    cx = 0.5; cy = 0.5; cz=0;
% 
%    [~,node2] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
% 	node2 = [node2;node2;node2];
%     dof2 = [1;2;3];
% 	val2 = [0;0;0];	
% 	bc_2 = [node2, dof2, val2];
% 	
%    cx = -0.5; cy = -0.5; cz=0;
%    
%    [~,node3] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
% 	node3 = [node3;node3;node3];
%     dof3 = [1;2;3];
% 	val3 = [0;0;0];	
% 	bc_3 = [node3, dof3, val3];
% 	
%    cx = 0.5; cy = -0.5; cz=0;
%    
%    [~,node4] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
% 	node4 = [node4;node4;node4];
%     dof4 = [1;2;3];
% 	val4 = [0;0;0];	
% 	bc_4 = [node4, dof4, val4];
% 	
%     bc.pDir = [ bc_1;
% 		        bc_2;
% 		        bc_3;
% 			    bc_4 ];

 %    1 apoyo fijo y 3 apoyos moviles

   cx = -0.5; cy = 0.5; cz=0;

    [~,node1] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
	node1 = [node1;node1;node1;node1];
    dof1 = [1;2;3;6];
	val1 = [0;0;0;0];	
	bc_1 = [node1, dof1, val1];
	
   cx = 0.5; cy = 0.5; cz=0;

   [~,node2] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
	node2 = [node2];
    dof2 = [3];
	val2 = [0];	
	bc_2 = [node2, dof2, val2];
	
   cx = -0.5; cy = -0.5; cz=0;
   
   [~,node3] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
	node3 = [node3];
    dof3 = [3];
	val3 = [0];	
	bc_3 = [node3, dof3, val3];
	
   cx = 0.5; cy = -0.5; cz=0;
   
   [~,node4] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
	node4 = [node4];
    dof4 = [3];
	val4 = [0];	
	bc_4 = [node4, dof4, val4];
	
    bc.pDir = [ bc_1;
		        bc_2;
		        bc_3;
			    bc_4 ];   

%% cond de neumman

   cx = 0; cy = 0; cz=0.2;
	
    [~,node] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);

    bc.pNeu = [node, 3, -150]; 
end