% clear all; close all; clc;
function [K,M,F] = solvper(mesh,a,c,f)
% %% PARAMETROS Y DATOS DEL PROBLEMA
%     % load problem datas
%     mesh = [];
%     params=[];
%     psi=[];
% 
% %       example= @perfilL; 
% %       example= @perfilcuadrado; 
% %       example=@hipercubo;
% %        example=@perfildobleT;
% %       example=@perfildobleT02;
%        example=@perfildobleT03;
% 
%     % load problem data
%     cd('examples')
%              [mesh, params, psi, bc, signatures, pdecoef, matprop] = example(mesh,psi,params);
%     cd ..

    % mesh and geometry parameters
    p = mesh.p;
    np = size(p,2);
    t = mesh.t;
    nt = size(t,2);
    area2=mesh.A;

nodos1 = t(1, :);
nodos2 = t(2, :);
nodos3 = t(3, :);
coord1 = p(:, nodos1);
coord2 = p(:, nodos2);
coord3 = p(:, nodos3);

lado13= coord1-coord3; 
lado23= coord2-coord3;
lado12= coord1-coord2;

% lado13=sum(lado13.^2);
% lado23=sum(lado23.^2);
% lado12=sum(lado12.^2);
% norma13=sqrt(lado13);
% norma23=sqrt(lado23);
% norma12=sqrt(lado12);


 % figure('Name','otro'); clf; set(1,'WindowStyle','docked');
 %        trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:)); 
 %        axis off


% determino el vector perpendicular

Czx= (lado13(2,:).*lado23(3,:)-lado13(3,:).*lado23(2,:)).^2;
Czy=(lado13(3,:).*lado23(1,:)-lado13(1,:).*lado23(3,:)).^2;
Czz= (lado13(1,:).*lado23(2,:)-lado13(2,:).*lado23(1,:)).^2;

norma_z= sqrt(Czx+Czy+Czz);

vers_z=[Czx./norma_z ; Czy./norma_z ; Czz./norma_z];

Cyx= (Czy.*lado13(3,:)-Czz.*lado13(2,:)).^2;
Cyy=(Czz.*lado13(1,:)-Czx.*lado13(3,:)).^2;
Cyz= (Czx.*lado13(1)-Czy.*lado13(1,:)).^2;

norma_y= sqrt(Cyx+Cyy+Cyz);

vers_y=[Cyx./norma_y ; Cyy./norma_y ; Cyz./norma_y];

Cxx= lado13(1,:).^2;
Cxy= lado13(2,:).^2;
Cxz= lado13(3,:).^2;

norma_x= sqrt(Cxx+Cxy+Cxz);

vers_x= [lado13(1,:)./norma_x; Cxy./norma_x; Cxz./norma_x];


x1l= (lado13(1,:).^2+lado13(2,:).^2+lado13(3,:).^2).^(0.5);
y1l= zeros(1,nt);
z1l= 0;

x2l= lado23(1,:).*vers_x(1,:)+lado13(2,:).*vers_x(2,:)+lado13(3,:).*vers_x(3);
y2l= lado23(1,:).*vers_y(1,:)+lado13(2,:).*vers_y(2,:)+lado13(3,:).*vers_y(3);
z2l= 0;

x3l= 0; 
y3l= 0;
z3l= 0;

x13= x1l;
x23= x2l;
y13= y1l;
y23= y2l;

area = norma_z./2;

area = reshape(area, 1, 1, nt);

x13 = reshape(x13, 1, 1, nt);
x23 = reshape(x23, 1, 1, nt);
y13 = reshape(y13, 1, 1, nt);
y23 = reshape(y23, 1, 1, nt);

B = zeros(2,3,nt);

B = [y23 -y13 (-y13-y23) ; -x23 x13 (x23-x13) ];

B= (1/(2*area(nt))).*B;

BT = permute(B, [2 1 3]);

ke= pagemtimes(BT, B);


mat= [2 1 1; 1 2 1 ; 1 1 2]/12;

mat = repmat(mat, [1 1 nt]);


mat= 1 * mat .* area;

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




% pchi = (psi<0); 
% tchi= pdeintrp(p, t, pchi);

Qe= [1;1;1];
Qe=repmat(Qe,[1 1 nt]);
f= reshape(f, 1, 1, nt);

re= (1/3).*area.*f.*Qe;

F1=zeros(np,1);
F2=zeros(np,1);
F3=zeros(np,1);

F1(nodos1,1)= reshape(re(1,1,:),[],1);
F2(nodos2,1)= reshape(re(2,1,:),[],1);
F3(nodos3,1)= reshape(re(3,1,:),[],1);

F= F1+F2+F3;

d= eigs(K,10,'smallestabs')

vector= ones(1,np);

vol= vector * M *vector'

end



     
     

    

  












