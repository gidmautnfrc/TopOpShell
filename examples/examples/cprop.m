function [pdecoef] = cprop(matprop)
%%% DEFINE PROPIEDADES CONSTITUTIVAS DEL PROBLEMA ENGLOBANDOLAS EN UN VECTOR
%%VALORES DE ENTRADA
% E= MODULO DE ELASTICIDAD
% nu = COEF DE POISSON
% h = ALTURA, PROPIEDAD GEOMETRICA
% gamma = NO ME ACUERDO TAMPOCO
%
% SALIDA
% Struct pdecoef.c_dkt 
%        pdecoef.c_andes
%
%%EL VECTOR C TIENE LA SIGUIENTE ESTRUCTURA
% c(1,1) = E; Modulo de elasticididad
% c(2,1) = Componente E1 del tensor constitutivo de DKT
% c(3,1) = Componente E2 del tensor constitutivo de DKT
%%
c=zeros(16,1);
%% Coeficientes de DKT

E=matprop.E0;
h=matprop.h0;
nu=matprop.nu0;


c(1,1) = E;
c(2,1) = (E.*h.^3)./(12.*(1-nu.^2));
c(3,1) = c(2,1)*nu;

%% Coeficientes de Andes

mu = matprop.mu_andes;
la = matprop.la_andes;
                        %%% Posiciones de c en pdecoef.ori
c(4,1) = 2*mu + la;       % 1 
c(5,1) = 0;               % 2 
c(6,1) = mu;              % 3 
c(7,1) = 0;               % 4 
c(8,1) = mu;              % 5 
c(9,1) = la;              % 6 
c(10,1) = 0;              % 7 
c(11,1) = mu;             % 8     
c(12,1) = 0;              % 9     
c(13,1) = 2*mu + la;      % 10  

pdecoef.c_dkt = c(2:3,:);
pdecoef.c_andes = c(4:13,1);


end