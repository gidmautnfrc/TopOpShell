function [mesh, params, psi0, bc, signatures, pdecoef, matprop] = memb(mesh,psi,params) %[mesh, pdecoef, matprop, params, bc, psi0]

    %% load g and b
    load mbb.mat;
    %% mesh generation
    
     if isempty(mesh)
        nsteps = 2;
        [p,e,t] = poimesh(g,20,10); t(4,:) = 1; 
        
        for i=1:2*nsteps
            [p,e,t] = refinemesh(g,p,e,t,'longest');            
        end
        
        mesh.p2d=p; %% Arreglo para que remalle con la opcion R
        psi0=-ones(length(unique(t(1:end-1,:))),1);
        
    else
            
        p = mesh.p2d; e = mesh.e; t = mesh.t; 
        [p,e,t,psi] = refinemesh(g,p,e,t,psi,'longest');                      
        [p,e,t,psi] = refinemesh(g,p,e,t,psi,'longest');  
        mesh.p2d=p;
        psi0 = psi;  
        
     end   
    
     %% Mesh struct

        [mesh]=meshstruct(mesh,p,e,t,g);
        
        %rotaciones de elementos

        [mesh]=rotation(mesh);
        
    
 
    %% material properties

     E = 210.0E9; 
%      E1 = 210.0E5; 
     nu = 0.3;  
     h = 1; 
     gamma = 1.0e-4; %% EN EL ORI LO CALCULA DENTRO DEL TOPDER, COMO K0/K1
                      % o en el solve lo calcula como E1/E0
%      matprop.E0 = E;
%      matprop.E1 = E1;
%      matprop.nu0 = nu;  

    matprop.E0 = E;
    matprop.gamma = gamma;
    matprop.nu0 = nu;
    matprop.h0 = h;
    matprop.E1 = gamma*E;
    
     %%Firma andes   

    [signatures]=andes_signature('OPT',nu);
    
    %% pde coeficients
    %DKT
    
    la = nu*E/((1+nu)*(1-2*nu)); 
    mu = E/(2*(1+nu)); % plane strain
    la = 2*mu*la/(la+2*mu); % plane stress
    
    matprop.la_dkt = la; 
    matprop.mu_dkt = mu;
    
    %ANDES
    la = nu*E/((1+nu)*(1-2*nu));
    mu = E/(2*(1+nu)); % plane strain
    la = 2*mu*la/(la+2*mu); % plane stress
    
    
%   matprop.la_andes = ((nu.*E.^2)./((1+nu).^2).*(1-2.*nu))./(((nu.*E)./(1+nu).*(1-2.*nu))+(E./(1+nu)));
    matprop.la_andes = la;
    matprop.mu_andes = E./(2.*(1+nu));
    
    %% Coeficiente c pde
    
    [pdecoef] = cprop(matprop);  
     
%     c = zeros(10,1);
%     c(1,:) = 2*mu + la; c(3,:) = mu;  c(5,:) = mu;
%     c(6,:) = la; c(8,:) = mu; c(10,:) = 2*mu + la;
% 
%     pdecoef.c = c;
%     pdecoef.a = zeros(4,1);
%     pdecoef.f = zeros(2,1);
%     pdecoef.b = b;        
    %% parameters
  
    % minimum allowed 'k'. Used by the line-search procedure
    params.kmin = 1.0E-3;
    % stop criterion
    params.stop = 1.0*pi/180;
    
    % method of penalization
    % (1) linear 
    % (2) exact quadratic
    % (3) augmented-lagrangian
    params.penalization = 1;
  
    % penalty parameter
    params.penalty = 10;     % method 1, 2 and 3
    params.volfrac = 0.3;   % method 2 and 3
    params.voleps  = 0.01;  % method 2 and 3
    params.auglag  = 0.5;   % method 3
    params.epsilon = 0.01;  % method 3       
        
    ghold  = 0;

    mesh.remesh = 'longest';
    mesh.ghold = ghold;    
    
  %% load-cases
    % Dirichlet conditions    
    bc.pDir = [];
    % bc.pDir = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal constraint 
  
    cx = 0;
    node = find(mesh.p(1,:) < cx + eps); node = node';
        
    npb = size(node,1); node = [node;node;node;node;node;node];
        
    dof = [ones(npb,1);2*ones(npb,1);3*ones(npb,1);4*ones(npb,1);5*ones(npb,1);6*ones(npb,1)];
        
    val = [zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1)];
    bc.pDir = [node, dof, val];
     
    
    % Neumman conditions
    bc.pNeu = [];  
    % lc.pNeu = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal load 

    cx = 40; cy = 10;
    node = intersect(find(p(1,:)==cx),find(p(2,:)==cy));
    bc.pNeu = [node, 2, -1.0];
     
    



    
   
  
end
