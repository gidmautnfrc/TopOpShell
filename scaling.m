function [scale,matprop,mesh,bc] = scaling (mesh,bc, matprop)
% SCALING
% Realiza un escalamientode de las variables a las al problema de placas 
%
% INPUT:
%   mesh:    Estructura de la malla original
%   matprop: Propiedades del material (debe incluir matprop.scaling_strategy)
%   bc:      Condiciones de borde originales
%
% OUTPUT:
%   mesh, matprop, bc: Estructuras modificadas según la estrategia ('A' o 'B')

%   Extraer dimensiones características del dominio
    % Define L_max como la máxima extensión en cualquier dirección
    Lx = max(mesh.p(1,:)) - min(mesh.p(1,:));
    Ly = max(mesh.p(2,:)) - min(mesh.p(2,:));
    Lz = max(mesh.p(3,:)) - min(mesh.p(3,:));
    L = [Lx Lz Lz];
    L= max(L);
    scale.L_c = L;

    mesh.p= mesh.p/scale.L_c;

 % % 2. Extraer cargas características de las condiciones de Neumann (bc.pNeu)
    F_max = max(abs(bc.pNeu));
    bc_m = F_max; % Fuerza característica de membrana 
    gc_b = F_max; % Fuerza característica de flexión 
    scale.b_c = bc_m;
    scale.g_c = gc_b;

    bc.pNeu = bc.pNeu/scale.b_c;
    
    % Parámetros físicos
    scale.E_c = matprop.E0;
    scale.gamma_c = matprop.gamma;
    scale.nu_c = matprop.nu0;
    scale.h_c = matprop.h0;
    scale.la_dkt_c = matprop.la_dkt ; 
    scale.mu_dkt_c = matprop.mu_dkt ;



    matprop.E0 = matprop.E0/scale.E_c;
    matprop.gamma =  matprop.gamma/scale.gamma_c;
    matprop.nu0 = matprop.nu0/scale.nu_c;
    matprop.h0 = matprop.h0/scale.h_c;
    matprop.la_dkt  = matprop.la_dkt/scale.la_dkt_c;
    matprop.mu_dkt =  matprop.mu_dkt/scale.mu_dkt_c;

    
end



   
    