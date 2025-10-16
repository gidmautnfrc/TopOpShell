
function [K,M,F] = assem_scalar_shell(mesh,c,a,f)

% extraccion de datos
p = mesh.p;
t = mesh.t;
np = size(p,2);
nt = size(t,2);

% extraccion los nodos por elemento de la matriz t
nodos1 = t(1, :);
nodos2 = t(2, :);
nodos3 = t(3, :);
coord1 = p(:, nodos1);
coord2 = p(:, nodos2);
coord3 = p(:, nodos3);

lado13= coord1-coord3; 
lado23= coord2-coord3;
                             
% determino versor x local
CXx= lado13(1,:);
CXy= lado13(2,:);
CXz= lado13(3,:);        
norma_CX= sqrt(CXx.^2+CXy.^2+CXz.^2);

vers_x= [lado13(1,:)./norma_CX; CXy./norma_CX; CXz./norma_CX];

% determino versor z local  CZ=lado13 x lado23 
CZx= lado13(2,:).*lado23(3,:)-lado13(3,:).*lado23(2,:);
CZy= lado13(3,:).*lado23(1,:)-lado13(1,:).*lado23(3,:);
CZz= lado13(1,:).*lado23(2,:)-lado13(2,:).*lado23(1,:);

norma_CZ= sqrt(CZx.^2+CZy.^2+CZz.^2);

vers_z=[CZx./norma_CZ ; CZy./norma_CZ ; CZz./norma_CZ];

%determino versor y local   CY= CZ x lado13
CYx= CZy.*lado13(3,:)-CZz.*lado13(2,:);
CYy= CZz.*lado13(1,:)-CZx.*lado13(3,:);
CYz= CZx.*lado13(1,:)-CZy.*lado13(1,:);

norma_CY= sqrt(CYx.^2+CYy.^2+CYz.^2);

vers_y=[CYx./norma_CY ; CYy./norma_CY ; CYz./norma_CY];

%coordenadas locales 

x13= norma_CX;
y13= zeros(1,nt);

x23= lado23(1,:).*vers_x(1,:)+lado13(2,:).*vers_x(2,:)+lado13(3,:).*vers_x(3);
y23= lado23(1,:).*vers_y(1,:)+lado13(2,:).*vers_y(2,:)+lado13(3,:).*vers_y(3);

area = norma_CZ./2;

area = reshape(area, 1, 1, nt);
x13 = reshape(x13, 1, 1, nt);
x23 = reshape(x23, 1, 1, nt);
y13 = reshape(y13, 1, 1, nt);
y23 = reshape(y23, 1, 1, nt);
f= reshape(f, 1, 1, nt);
B = zeros(2,3,nt);

B = [y23 -y13 (y13-y23) ; -x23 x13 (x23-x13) ];

B= (1./(2*area)).*B;

BT = permute(B, [2 1 3]);

ke=c.*pagemtimes(BT, B).*area;

mat= [2 1 1; 1 2 1 ; 1 1 2]/12;

mat = repmat(mat, [1 1 nt]);

mat= mat .* (area.*a) ;

M = sparse(np,np);

for i = 1:3
    for j = 1:3
        M = M + sparse(t(i,:),t(j,:),reshape(mat(i,j,:),1,[]),np,np);
    end
end

K = sparse(np,np);
for i = 1:3
    for j = 1:3
        K = K + sparse(t(i,:),t(j,:),reshape(ke(i,j,:),1,[]),np,np);
    end
end

Q=repmat([1;1;1],[1 1 nt]);

re= (1/3).*area.*f.*Q;

F = sparse(np, 1);
F = F + sparse(t(1, :), 1, reshape(re(1, 1, :), 1, []), np, 1);
F = F + sparse(t(2, :), 1, reshape(re(2, 1, :), 1, []), np, 1);
F = F + sparse(t(3, :), 1, reshape(re(3, 1, :), 1, []), np, 1);
F = full(F);
end



     
     

    

  












