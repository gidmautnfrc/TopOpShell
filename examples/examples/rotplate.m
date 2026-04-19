function [mesh, params, psi0, bc, signatures, pdecoef, matprop] = rotplate(mesh,psi,params)
    %%cargar malla
   
    run('platerot.m') % Coordenadas y malla = [m]
    t=msh.TRIANGLES';

    if isempty(mesh)
               psi0=-ones(length(unique(t(1:end-1,:))),1);   
    else
        psi0 = psi;          
    end   
    
    %% Mesh struct
        
        [mesh]=meshstruct2(msh);
    
    %rotaciones de elementos

        [mesh]=rotation(mesh);
    
    %% Material properties
  
    E = 70.0E+9; 
    gamma = 1.0e-4;
    nu = 0.3;
    h = 1; 
    
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
        params.penalty = 15; 
        
    %% boundary conditions
    
    bc.pNeu = []; bc.pDir = []; 
    x = mesh.p(1,:); y = mesh.p(2,:); z = mesh.p(3,:);
    eps = 0.01*sqrt(min(mesh.A));

        cz = 0.0;     
        node = find(z <  eps); node = node';
        npb = size(node,1); node = [node;node;node;node;node;node];
        dof = [ones(npb,1);2*ones(npb,1);3*ones(npb,1);4*ones(npb,1);5*ones(npb,1);6*ones(npb,1)];
        val = [zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1)];
        bc.pDir = [node, dof, val];
        

%% cond de neumman
    cx = 0.5; cy = 1; cz=0.5;
    node = min(find((x-cx).^2 + (y-cy).^2 <= eps^2)); 
    % 1  --> Fx
    % 2  --> Fy
    % 3  --> Fz
    % 4  --> Mx
    % 5  --> My
    % 6  --> Mz
    bc.pNeu = [node, 1,-100;
               node, 3,-100]; 
end
