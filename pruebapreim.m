clear all;close all;clc;

%% PARAMETROS Y DATOS DEL PROBLEMA
    % load problem datas
    mesh = [];
    params=[];
    psi=[];
    
%       example= @perfilL; 
%       example= @perfilcuadrado; 
%       example=@hipercubo;
%        example=@perfildobleT;
%       example=@perfildobleT02;
       example=@perfildobleT03;

    % load problem data
    cd('examples')
             [mesh, params, psi, bc, signatures, pdecoef, matprop] = example(mesh,psi,params);
    cd ..

p= mesh.p; t=mesh.t;
% alpha= params.alpha * per_defined;
alpha=1;
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

centro = [0.5; 0.25; 0.25]; 
r= 0.1; 

restriccion = ((p(1,:) - 0.5).^2 + (p(2,:) - 0.25).^2 + (p(3,:) - 0.25).^2) <= r^2;
psi(restriccion)=1;

pchi = (psi<0); %funcion caracteristica 



tchi= pdeintrp(p, t, pchi); %funcion caracteristica interpolada en el baricentro de cada elemento       
ep =2.1*h;
c = ep^2;
a = 1;
f = tchi;
m=0;


 % topplot=1*(psi>0);
 %    tchilogical=find(tchi>0);
 %    tplot= t(:,tchilogical);
 %    tchiplot=-tchi(tchilogical);    
    
    % figure(3); clf; set(3,'WindowStyle','docked');
    % trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:), tchi,'LineStyle','none');
    % colormap
    % axis equal

[K,M,F] = assem_scalar_shell(mesh,c,a,f);

% [K,M,F] = assema(p,t,c,a,f);
K=K+M;
vNew = K \ F;


    figure(3); clf; set(3,'WindowStyle','docked');
    trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:), vNew,'LineStyle','none');
    axis equal
 
    figure(4); clf; set(3,'WindowStyle','docked');
    trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:), tchi,'LineStyle','none');
    axis equal
    colorbar

perimeter=(2/ep)*(1-vNew)'*M*pchi
% perimeter=(perimeter-8.979114584668373e+01)/2


circulo= 2*pi*r

% dtper= (1/ep)*(1-2*vNew);
% 
% dtper= dtper/Per0 * alpha;

% end