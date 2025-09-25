clear all; close all; clc;
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

    % mesh and geometry parameters
    
    p = mesh.p;
    np = size(p,2);
%     e = mesh.e;
    t = mesh.t;
    nt = size(t,2);

nodos1 = t(1, :);
nodos2 = t(2, :);
nodos3 = t(3, :);
coord1 = p(:, nodos1);
coord2 = p(:, nodos2);
coord3 = p(:, nodos3);
lado12 = coord1 - coord2;
lado13 = coord1 - coord3;
% Extraer los componentes de los vectores de los lados
x12 = lado12(1, :);
y12 = lado12(2, :);
z12 = lado12(3, :);

x13 = lado13(1, :);
y13 = lado13(2, :);
z13 = lado13(3, :);

% Calcular las componentes del vector de producto cruz para cada elemento
% Componente x del resultado:
Cx = y12 .* z13 - z12 .* y13;

% Componente y del resultado:
Cy = z12 .* x13 - x12 .* z13;

% Componente z del resultado:
Cz = x12 .* y13 - y12 .* x13;

% Unir los componentes para formar la matriz de vectores de producto cruz
% El resultado 'cross_prod' es una matriz 3xN donde cada columna es el producto cruz
cross_prod = [Cx; Cy; Cz];
squared_elements = cross_prod.^2;
sum_of_squares = sum(squared_elements, 1);
modulos = sqrt(sum_of_squares);
area = modulos./2;

area = reshape(area, 1, 1, nt);

mat= [2 1 1; 1 2 1 ; 1 1 2];
mat = repmat(mat, [1 1 nt]);


mat= 0.003 * mat .* area;
M = sparse(np,np);
for i = 1:3
    for j = 1:3
        M = M + sparse(t(i,:),t(j,:),reshape(mat(i,j,:),1,[]),np,np);
    end
end



    % nodo1=t(1,:);
    % nodo2=t(2,:);
    % nodo3=t(3,:);
% 
% i = []; % Índices de Fila Global
% j = []; % Índices de Columna Global
% v = []; % Valores
% 
% for e=size(area)
% 
%     M_e= area(e)*mat;
%     nodo1=t(1,e);
%     nodo2=t(2,e);
%     nodo2=t(3,e);
% 
%     i_e= (nodo1,nodo1,nodo1,nodo2,nodo2,nodo2,nodo3,nodo3,nodo3);
%     j_e= (nodo1,nodo2,nodo3,nodo1,nodo2,nodo3,nodo1,nodo2,nodo3);
%     v_e= (M_e(1,1), M_e(1,2), M_e(1,3), M_e(2,1), M_e(2,2), M_e(2,3), M_e(3,1), M_e(3,1),M_e(3,3));
% 
%     ia(:,e)= (i_e)';
%     ja(:,e)= (j_e)';
%     va(:,e)= (v_e)';
% end
% 
% 
% i= ia(:)';
% j= ja(:)';
% v= va(:)';
% M = sparse(i,j,v,size(p,1),size(p,1))

     
     

    

  












