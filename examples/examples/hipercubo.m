function [mesh, params, psi0, bc, signatures, pdecoef, matprop] = hipercubo(mesh,psi,params)
    %%cargar malla
	
	mallas = {'hipercubo_mallag.mphtxt','hipercubo_mallam.mphtxt','hipercubo_mallaf.mphtxt','hipercubo_mallaf02.mphtxt'}; % 'hipercubomin.mphtxt','hipercubo.mphtxt'
	
     if isempty(mesh);
	 
	        idx = 1;
            mesh = getComsolmesh(mallas{idx});
               psi0 = -ones(length(unique(mesh.t(1:end-1,:))),1);   
    else
	
		idx = mesh.index;
		idx = idx + 1;
		
		p_old = mesh.p';
		
		mesh = getComsolmesh(mallas{idx});
		
		p_new = mesh.p';
        
        I = scatteredInterpolant(p_old, psi);
        psi0 = I(p_new);          
    end   
    
	mesh.index = idx;
	
    %% Mesh struct
        
        [mesh] = meshstruct3(mesh);
    
	
    %rotaciones de elementos

        [mesh] = rotation(mesh);
    
    %% Material properties
  
    E = 70.0E+9;  %% MPa
    gamma = 1.0e-4;
    nu = 0.3; 
    h = 0.003; %% m
    
    matprop.E0 = E;
    matprop.gamma = gamma;
    matprop.nu0 = nu;
    matprop.h0 = h;
    matprop.E1 = gamma*E;
    
    
    %DKT
    
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
        params.penalty = 2000; 
        %params.penalty = 50000;
    %% boundary conditions
    
    bc.pNeu = []; bc.pDir = []; 
    x = mesh.p(1,:); y = mesh.p(2,:); z = mesh.p(3,:);
    eps = 0.01*sqrt(min(mesh.A));

        cx = 0.10;     
        node = find(x > cx - eps); node = node';
        npb = size(node,1); node = [node;node;node;node;node;node];
        dof = [ones(npb,1);2*ones(npb,1);3*ones(npb,1);4*ones(npb,1);5*ones(npb,1);6*ones(npb,1)];
        val = [zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1)];
        bc.pDir = [node, dof, val];

%% cond de neumman
    cx = 0.0; cy = 0.025; cz=0.025;
	
	eps1 = 1.01 * min((y - cy).^2 + (z - cz).^2);
    node1 = min(find((y-cy).^2 + (z-cz).^2 <= eps1)); 
%     cx = 0.0; cy = 0.0; cz=0.0;
%     node2 = min(find((x-cx).^2 + (y-cy).^2 <= eps^2)); 

    bc.pNeu = [node1, 3,-150]; 
end