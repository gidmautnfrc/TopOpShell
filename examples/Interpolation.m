function [psi] = Interpolation(psi,mesh,msh)

t_old = mesh.t; % tiene un elemento por columna
p_old = mesh.p; % Tiene un nodo por columna

%Armo vectores con los indices de las columnas de la matriz de conectividad
id1 = t_old(1,:); id2 = t_old(2,:); id3 = t_old(3,:);

x1 = p_old(1,id1); x2 = p_old(1,id2); x3 = p_old(1,id3);
y1 = p_old(2,id1); y2 = p_old(2,id2); y3 = p_old(2,id3);
z1 = p_old(3,id1); z2 = p_old(3,id2); z3 = p_old(3,id3);

x_g = (x1+x2+x3)/3; y_g = (y1+y2+y3)/3; z_g = (z1+z2+z3)/3;

% Hasta aca tome la matriz de conectividad vieja y las coordenadas
% y calcule cordenadas baricentricas

%% Ahora lo que necesito son las nuevas coordenadas porque tengo
% que interpolar cada una de esas cordenadas nuevas.
% Esto lo consigo tomando los valores de la nueva malla msh

p_new = msh.POS'; % Tiene un nodo por columna
% t_new = msh.TRIANGLES'; % tiene un elemento por columna

num = size(p_new,2); %p es de la nueva malla

psi_new = zeros(num, 1);

for i = 1:num

x_new = p_new(1,i); y_new = p_new(2,i); z_new = p_new(3,i);

[~,id] = min((x_new - x_g).^2 + (y_new - y_g).^2 + (z_new - z_g).^2); %Me devuelve el indice del elemento cuya distancia es la menor

% Ya encontre el elemento que contiene al punto
% ahora necesito calcular las funciones de forma

v1 = p_old(:,t_old(1,id)); v2 = p_old(:,t_old(2,id)); v3 = p_old(:,t_old(3,id));
P = [x_new; y_new; z_new];

A = 0.5*norm(cross(v2-v1,v3-v1));
A1 = 0.5*norm(cross(P-v2,v3-v2));
A2 = 0.5*norm(cross(P-v1,v3-v1));
A3 = 0.5*norm(cross(P-v1,v2-v1));

%Armo el vector de funciones de forma
N = [A1/A; A2/A; A3/A];

psi_v = psi([t_old(1,id); t_old(2,id); t_old(3,id)]);

psi_new(i) = dot(N,psi_v); %Una interpolacion lineal tambien se puede expresar como un producto escalar

end

psi = psi_new;

end