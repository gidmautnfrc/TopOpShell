function [dtper,perimeter] = TopDerPer(mesh,psi,params,per_defined,Per0)
 
p= mesh.p; t=mesh.t;
alpha= params.alpha * per_defined;

%% Determain  mean largest length 
% 1. Obtener los índices de los nodos para cada vértice de los triángulos
%    Se leen las filas de 't'
nodos1 = t(1, :);
nodos2 = t(2, :);
nodos3 = t(3, :);

% 2. Obtener las coordenadas (x,y,z) de cada vértice
%    Se leen las columnas de 'p' correspondientes a los índices
coord1 = p(:, nodos1);
coord2 = p(:, nodos2);
coord3 = p(:, nodos3);

% 3. Calcular los vectores que definen los lados de cada triángulo
lado12 = coord1 - coord2;
lado23 = coord2 - coord3;
lado31 = coord3 - coord1;

% 4. Calcular el cuadrado de las longitudes de los lados
%    Se suma a lo largo de la dimensión 1 (filas: dx^2 + dy^2 + dz^2)
longitudes_sq_12 = sum(lado12.^2, 1);
longitudes_sq_23 = sum(lado23.^2, 1);
longitudes_sq_31 = sum(lado31.^2, 1);

% 5. Encontrar la longitud máxima al cuadrado para cada triángulo
%    Se apilan los vectores de longitud y se busca el máximo en la dimensión 1 (columna por columna)
longitudes_max_sq = max([longitudes_sq_12; longitudes_sq_23; longitudes_sq_31], [], 1);

% 6. Calcular la raíz cuadrada para obtener la longitud máxima real
l_max = sqrt(longitudes_max_sq);
longitud_media_maxima = mean(l_max);

% Determinacion de h

h= longitud_media_maxima;

%resolucion de  vNeu

h = mean(l_max);
pchi = (psi<0); %funcion caracteristica 
tchi= pdeintrp(p, t, pchi); %funcion caracteristica interpolada en el baricentro de cada elemento       
ep =4*h;
c = ep^2;
a = 1;
f = tchi;
m=0;

[K,M,F] = solvper(mesh,a,c,f);






% [K,M,F] = assema(p,t,c,a,f);
K=K+M;
vNew = K \ F;

%perimeter = (1/ep)*v'*M*pchi;
perimeter=(2/ep)*(1-vNew)'*M*pchi;

dtper= (1/ep)*(1-2*vNew);

dtper= dtper/Per0 * alpha;

end