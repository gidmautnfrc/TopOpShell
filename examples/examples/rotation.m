function [mesh]=rotation(mesh);
%**************************************************************************
%Traslacion y rotacion de matrices de un sistema global a un sistema local%
%**************************************************************************
%
% DESCRIPTION
%
%  Calcula desde los puntos en un sistema de coordenadas global, la
%  traslacion y rotacion a un sistema de coordenadas local, tomando como
%  origen el primer punto del conjunto de puntos realacionados al sistema
%  global. Realizando una transformacion (x,y,z) a (x',y').
%
% INFORMACION
%
%  Variables de entrada: mesh.p = Puntos del plano en el sistema de coordenadas
%                            global
%                        mesh.t = Conectividad
%   
%  Salida: p = Puntos del plano en el sistema de coordenadas global (x,y,z)
%
%          plnodos (1,2,3) = Puntos de los nodos de los elementos en el
%                              sistema de coordenadas local (x,y)
%
%          lamdagtol = Matriz de cosenos directores
%
% HISTORY
% E.D. Vigna     08/2019: Code implementation.
% E.D. Vigna     19/2022: Modify entry and out for struct mesh
%
%**************************************************************************
% %
% HIPOTESIS
%
% Se consideran tres puntos de p, los que conforman cada elemento, en
% sentido antihorario. Estos puntos se deben definir en el sistema de
% coordenadas globales.
% En el sistema coordenado local el punto 1, corresponde al origen.
% Calculando la traslacion del sistema global al local referenciando dicho
% punto y luego realizando la rotacion.
%
%                .3
%               /|
%              / |
%             /  |
%            /   |
%           /    |
%          1.----2.
%
%
%
%**************************************************************************

%% ROTACION
% Obtener valores del struct mesh 
    p=mesh.p;
    t=mesh.t;
% LAMDAS PARA CADA ELEMENTO

uvec =  p(:,t(2,:)) - p(:,t(1,:));
vvec =  p(:,t(3,:)) - p(:,t(1,:));

% VERSORES COORDENADAS LOCALES

eigenjlocalaux = vvec ./ (((vvec(1,:).^2+(vvec(2,:)).^2+(vvec(3,:).^2)).^(1/2)));
eigenilocal = uvec ./ (((uvec(1,:).^2+(uvec(2,:)).^2+(uvec(3,:).^2)).^(1/2)));

eigenklocal = [(eigenilocal(2,:).*eigenjlocalaux(3,:)-eigenilocal(3,:).*eigenjlocalaux(2,:)) ;...
    (eigenilocal(3,:).*eigenjlocalaux(1,:)-eigenilocal(1,:).*eigenjlocalaux(3,:)) ;...
    (eigenilocal(1,:).*eigenjlocalaux(2,:)-eigenilocal(2,:).*eigenjlocalaux(1,:)) ];

eigenklocal = eigenklocal./(((eigenklocal(1,:)).^2+(eigenklocal(2,:)).^2+(eigenklocal(3,:)).^2).^(1/2));

eigenjlocal = [(eigenklocal(2,:).*eigenilocal(3,:)-eigenklocal(3,:).*eigenilocal(2,:)) ;...
    (eigenklocal(3,:).*eigenilocal(1,:)-eigenklocal(1,:).*eigenilocal(3,:)) ;...
    (eigenklocal(1,:).*eigenilocal(2,:)-eigenklocal(2,:).*eigenilocal(1,:)) ];

% COSENOS DIRECTORES SISTEMAS DE COORDENADAS LOCAL GLOBAL

lambdaxlxg = eigenilocal(1,:);
lambdaxlyg = eigenilocal(2,:);
lambdaxlzg = eigenilocal(3,:);

lambdaylxg = eigenjlocal(1,:);
lambdaylyg = eigenjlocal(2,:);
lambdaylzg = eigenjlocal(3,:);

lambdazlxg = eigenklocal(1,:);
lambdazlyg = eigenklocal(2,:);
lambdazlzg = eigenklocal(3,:);


%% TRASLACION Y OBTENCION DE COORDENADAS LOCALES DE P

% TRASLACION NODOS

ptrasnod1 = p(:,t(1,:)) - p(:,t(1,:));
ptrasnod2 = p(:,t(2,:)) - p(:,t(1,:));
ptrasnod3 = p(:,t(3,:)) - p(:,t(1,:));

% ROTACION NODOS

xlnod1 = ptrasnod1(1,:).*lambdaxlxg +  ptrasnod1(2,:).*lambdaxlyg +  ptrasnod1(3,:).*lambdaxlzg;
ylnod1 = ptrasnod1(1,:).*lambdaylxg +  ptrasnod1(2,:).*lambdaylyg +  ptrasnod1(3,:).*lambdaylzg;
zlnod1 = ptrasnod1(1,:).*lambdazlxg +  ptrasnod1(2,:).*lambdazlyg +  ptrasnod1(3,:).*lambdazlzg;

xlnod2 = ptrasnod2(1,:).*lambdaxlxg +  ptrasnod2(2,:).*lambdaxlyg +  ptrasnod2(3,:).*lambdaxlzg;
ylnod2 = ptrasnod2(1,:).*lambdaylxg +  ptrasnod2(2,:).*lambdaylyg +  ptrasnod2(3,:).*lambdaylzg;
zlnod2 = ptrasnod2(1,:).*lambdazlxg +  ptrasnod2(2,:).*lambdazlyg +  ptrasnod2(3,:).*lambdazlzg;

xlnod3 = ptrasnod3(1,:).*lambdaxlxg +  ptrasnod3(2,:).*lambdaxlyg +  ptrasnod3(3,:).*lambdaxlzg;
ylnod3 = ptrasnod3(1,:).*lambdaylxg +  ptrasnod3(2,:).*lambdaylyg +  ptrasnod3(3,:).*lambdaylzg;
zlnod3 = ptrasnod3(1,:).*lambdazlxg +  ptrasnod3(2,:).*lambdazlyg +  ptrasnod3(3,:).*lambdazlzg;

%AGRUPACION DE DATOS

plnod1 = [ xlnod1 ; ylnod1; zlnod1];
plnod2 = [ xlnod2 ; ylnod2; zlnod2];
plnod3 = [ xlnod3 ; ylnod3; zlnod3];

%COORDENADAS PARAMETRICAS LOCALES

x1 = [plnod1(1,:)];
x2 = [plnod2(1,:)];
x3 = [plnod3(1,:)];

y1 = [plnod1(2,:)];
y2 = [plnod2(2,:)];
y3 = [plnod3(2,:)];

%% AREA

A = ((((uvec(2,:).*vvec(3,:))-(uvec(3,:).*vvec(2,:))).^2+((uvec(3,:).*vvec(1,:))-(uvec(1,:).*vvec(3,:))).^2+((uvec(1,:).*vvec(2,:))-(uvec(2,:).*vvec(1,:))).^2).^(1/2)).*(1/2);

%% LONG

long12 = ((x2-x1).^2+(y2-y1).^2).^0.5;
long23 = ((x3-x2).^2+(y3-y2).^2).^0.5;
long31 = ((x3-x1).^2+(y3-y1).^2).^0.5;

% Sin raices para ANDES

l12 = ((x2-x1).^2+(y2-y1).^2);
l23 = ((x3-x2).^2+(y3-y2).^2);
l31 = ((x3-x1).^2+(y3-y1).^2); 

%% OUT struct mesh

mesh.lambda.lambdaxgxl = lambdaxlxg;
mesh.lambda.lambdaxgyl = lambdaxlyg;
mesh.lambda.lambdaxgzl = lambdaxlzg;

mesh.lambda.lambdaygxl = lambdaylxg;
mesh.lambda.lambdaygyl = lambdaylyg;
mesh.lambda.lambdaygzl = lambdaylzg;

mesh.lambda.lambdazgxl = lambdazlxg;
mesh.lambda.lambdazgyl = lambdazlyg;
mesh.lambda.lambdazgzl = lambdazlzg;

%Local to global

mesh.lambda.lambdaxlxg = lambdaxlxg;
mesh.lambda.lambdaxlyg = lambdaxlyg;
mesh.lambda.lambdaxlzg = lambdaxlzg;
mesh.lambda.lambdaylxg = lambdaylxg;
mesh.lambda.lambdaylyg = lambdaylyg;
mesh.lambda.lambdaylzg = lambdaylzg;
mesh.lambda.lambdazlxg = lambdazlxg;
mesh.lambda.lambdazlyg = lambdazlyg;
mesh.lambda.lambdazlzg = lambdazlzg;

%Local Coordinates

mesh.localcoordinates.x1 = x1;
mesh.localcoordinates.x2 = x2;
mesh.localcoordinates.x3 = x3;
mesh.localcoordinates.y1 = y1;
mesh.localcoordinates.y2 = y2;
mesh.localcoordinates.y3 = y3;

%Long without root for andes

mesh.long.l12 = l12;
mesh.long.l23 = l23;
mesh.long.l31 = l31;

%Long side elements

mesh.long.long12 = long12;
mesh.long.long23 = long23;
mesh.long.long31 = long31;   

%area
mesh.A=A;
end       