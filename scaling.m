function [mesh, params, psi0, bc, signatures, pdecoef, matprop] = scaling (scale)
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
    L = max(Lx, Ly, Lz); 
    scale.L_c = L;

    p= p/scale.L_c;

 % % 2. Extraer cargas características de las condiciones de Neumann (bc.pNeu)
    F_max = max(abs(bc.pNeu(:,3)));
    bc_m = F_max; % Fuerza característica de membrana 
    gc_b = F_max; % Fuerza característica de flexión 
    scale.m_c = bc_m;
    scale.b_c = gc_b;

    bc.pNeu = bc.pNeu/m_c;
    
    % Parámetros físicos
    scale.E_c = matprop.E0;
    scale.gamma_c = matprop.gamma;
    scale.nu_c = matprop.nu0;
    scale.h_c = matprop.h0;
    scale.la_dkt_c = matprop.la; 
    scale.mu_dkt_c = matprop.mu;


    matprop.E0 = matprop.E0/scale.E_c;
    matprop.gamma =  matprop.gamma/scale.gamma_c;
    matprop.nu0 = matprop.nu0/scale.nu_c;
    matprop.h0 = matprop.h0/scale.h_c;
    matprop.la = matprop.la/scale.la_dkt_c;
    matprop.mu =  matprop.mu/scale.mu_dkt_c;

    
end



   
    