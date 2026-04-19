function [mesh, params, psi0, bc, signatures, pdecoef, matprop] = scaling ()
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
    scale.L_max = L;

 % % 2. Extraer cargas características de las condiciones de Neumann (bc.pNeu)
    % Asumiendo que la 3ra columna de bc.pNeu contiene los valores de fuerza
    % Se requiere una lógica robusta para separar membrana (dof 1,2) y flexión (dof 3)
    % Aquí se esquematiza el concepto:
    F_max = max(abs(bc.pNeu(:,3))); 
    
    % ATENCIÓN: Debes programar la extracción correcta de b_c y g_c basada en 
    % los grados de libertad asociados a la membrana y a la flexión.
    bc_m = F_max; % Fuerza característica de membrana (Placeholder)
    gc_b = F_max; % Fuerza característica de flexión (Placeholder)
    
    matprop.force_m = bc_m;
    matprop.force_b = gc_b;
    
    % Parámetros físicos
    h = matprop.h0;
    mu = matprop.mu_dkt; % Módulo de corte
    
    % 3. Bifurcación según la estrategia de escalamiento
    switch upper(matprop.scaling_strategy)
        case 'A'
            % ESTRATEGIA A: Adimensionalización fuerte (Forma débil)
            
            % Escalar geometría
            mesh.p = mesh.p / L;
            
            % Escalar cargas de contorno (Adimensionalización de F)
            % Pseudo-código: bc.pNeu(es_membrana, 3) = bc.pNeu(es_membrana, 3) / bc_m;
            % Pseudo-código: bc.pNeu(es_flexion, 3) = bc.pNeu(es_flexion, 3) / gc_b;
            
            % Computar factores de escala de desplazamiento
            matprop.wc = (12 * L^4 * gc_b) / (h^3 * mu);
            matprop.uc = matprop.wc * (gc_b / bc_m);
            
            % Notificamos que Pi2 no se usa en la estrategia A
            matprop.Pi2 = 1.0; 
            
        case 'B'
            % ESTRATEGIA B: Escalamiento del funcional de costo (Penalización)
            
            % La malla y las cargas se mantienen con dimensiones físicas.
            matprop.uc = 1.0;
            matprop.wc = 1.0;
            
            % Computar el número adimensional Pi_2
            matprop.Pi2 = (gc_b * L^2) / (bc_m * h^2);
            
        case 'NONE'
            % Sin escalamiento (Comportamiento original con mal condicionamiento)
            matprop.uc = 1.0;
            matprop.wc = 1.0;
            matprop.Pi2 = 1.0;
            matprop.L_char = 1.0;
            
        otherwise
            error('Estrategia de escalamiento no reconocida. Use ''A'', ''B'', o ''NONE''.');
    end
end



   
    