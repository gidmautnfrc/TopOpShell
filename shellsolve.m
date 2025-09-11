function [U,F] = shellsolve(mesh, pdecoef, matprop, signatures, bc,psi)
% mesh parameters
    p = mesh.p; t = mesh.t; np = size(p,2);
% material properties
    gamma = matprop.gamma; % contrast  
    tgamma = pdeintrp(p,t,(psi<0) + gamma*(psi>=0)); 

% pde coefficients    
   c_dkt0 = pdecoef.c_dkt;
   c_andes0= pdecoef.c_andes;
   
   pdecoef.c_dkt = c_dkt0*tgamma;
   pdecoef.c_andes= c_andes0*tgamma;
%%Matriz de rigidez
    
    [K_SHELL]= assemashell(mesh,pdecoef, matprop, signatures);
    
     K = K_SHELL;
    
%% Matriz F

    F=sparse(mesh.np*mesh.N,1);
    
    %% condiciones de borde
    
        % apply Dirichlet boundary condition
        np=mesh.np;
        eq = bc.pDir(:,1) + (bc.pDir(:,2)-1)*np; neq = size(eq,1);
        %eq = nodo + gdl-1*numero de nodos
        K(eq,:) = 0.0; K(:,eq) = 0.0; K(eq,eq) = speye(neq);

        % apply Neumann boundary condition
        eq = bc.pNeu(:,1) + (bc.pNeu(:,2)-1)*np; 
        F(eq) = F(eq) + bc.pNeu(:,3);   

    
    %% Calculo desplazamientos
    
    U_shell = K \ F;   
    [U] = disp_filter(mesh, U_shell);
    
end