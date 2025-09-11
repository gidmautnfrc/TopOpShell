%**************************************************************************
% Topological Derivative for the Compliance  
%**************************************************************************
% DESCRIPTION
% Computes the topological derivative
% 
% INPUT
% mesh:    pdetool mesh struct
% U:       pdetool solution
% matprop: material properties struct
% psi:     level-set function
% params:  topology optimization parameters struct
% pdecoef: pdetool coeficiets struct
%
% OUTPUT
% dtE: compliance T. D. -> dtE = -DT (bulk) and dt = DT (inclusion)
%
% HISTORY
% A.A. Novotny     05/2012: code implementation.
% C.G. Lopes       07/2017: final code updating.
% R.B. dos Santos  07/2017: final code updating.
% S.M. Giusti      08/2018: implementation of RM plate.
%**************************************************************************

function [dt] = tdshell(mesh,U,pdecoef, matprop,signatures,psi)
%% Preliminaries generales
    np=mesh.np;
    p=mesh.p;
    t=mesh.t;
    h=matprop.h0;
    nu=matprop.nu0;
    la=matprop.la_andes;
    mu=matprop.mu_andes;

%% Problema de membrana
        % Preliminaries
            
    gamma = matprop.gamma;
    tgamma = pdeintrp(p,t,(psi<0)+gamma*(psi>=0)); %%implementar interpolacion   
    id=[1;1;0];
    
        alpha = (la+mu)/mu;
        beta = (la+3*mu)/(la+mu);
        e = e_andes(mesh,U,h,signatures); % strain     
        e(3,:)=1/2.* e(3,:); 
        divu = e(1,:) + e(2,:);  % e(1,:) = exx ; e(2,:) = eyy ; e(3,:) = exy de ANDES
        s = la*id*divu + 2*mu*e; % Nominal stress                
        tgamma3 = [tgamma;tgamma;tgamma]; s = s.*tgamma3; % Effective stress  
        
        %   Incluiding the body force contribution
        %         f = f0;
        %         f=0;
        %         tU = pdeintrp(p,t,U); 
        %         dtf = (f(1,:).*tU(1,:)+f(2,:).*tU(2,:));        

        % Topological derivative at the bulk phase
        coef0 = 0.5*((1.0-gamma)/(1+gamma*beta)); 
        coef1 = coef0*(1+beta); 
        coefx = (1.0-gamma)/(1.0+gamma*alpha);         
        coef2 = 0.5*coef0*(alpha-beta)*coefx;         

        %             Acomodado como el arreglo de matlab, donde e(1)=exx ;
        %             e2=exy; e3=eyy; entoontes donde e(1,:) = eandes(1,:)
        %                                             e(2,:) = eandes(3,:)
        %                                             e(3,:) = eandes(2,:)

           dte_m = coef1*((s(1,:).*e(1,:)+2*s(3,:).*e(3,:)+s(2,:).*e(2,:))) ...
            + coef2*((s(1,:)+s(2,:)).*(e(1,:)+e(2,:)));% + dtfe ;  
                                        
        % Topological derivative at the inclusion

        gamma = 1.0/gamma; 
        coef0 = 0.5*((1-gamma)/(1+gamma*beta)); coef1 = coef0*(1+beta);
        coefx = (1.0-gamma)/(1.0+gamma*alpha);
        coef2 = 0.5*coef0*(alpha-beta)*coefx;

        dti_m = coef1*((s(1,:).*e(1,:)+2*s(3,:).*e(3,:)+s(2,:).*e(2,:))) ...
            + coef2*((s(1,:)+s(2,:)).*(e(1,:)+e(2,:)));% + dtfi; 
        
%% Problema de placas
% Preliminaries
    
gamma = matprop.gamma;
tgamma = pdeintrp(p,t,(psi<0)+gamma*(psi>=0));
    
h0 = matprop.h0;
h3 = h0*h0*h0;
la=matprop.la_dkt;
mu=matprop.mu_dkt;

alpha = (la+mu)/mu;
beta = mu/(3.0*mu+2.0*la);
nu = la/(2*mu+la);
E = 2*mu*(1+nu);
tE = E*tgamma;
s = getstress_plate(mesh,U,matprop); % nominal stress      CAMBIAR A PLATES NO COMP
tgamma3 = [tgamma;tgamma;tgamma]; 
s = s.*tgamma3; % effective stress
aa = s(1,:).*s(1,:) + 2*s(2,:).*s(2,:) + s(3,:).*s(3,:); % s.s
bb = (s(1,:) + s(3,:)).*(s(1,:) + s(3,:)); % (tr(s))^2

% Topological derivative at the bulk phase
% coef0 = 6.0*(1.0 - gamma)./(tE*h3);  %% EN EL ORI EL COEF0 ES 12*1-GAMMA...
 coef0 = 12.0*(1.0-gamma)./(tE*h3); 
coef1 = 4.0*(alpha*beta)/(1.0 + beta*gamma);
coef2 = 1.0/(1.0 + alpha*gamma) - 2.0*(alpha*beta)/(1.0 + beta*gamma);
dte_b = coef0.*(coef1*aa + coef2*bb);

% Topological derivative at the inclusion
gamma = 1.0/gamma;
% coef0 = 6.0*(1.0 - gamma)./(tE*h3);%% EN EL ORI EL COEF0 ES 12*1-GAMMA...
coef0 = 12.0*(1.0-gamma)./(tE*h3); %% AcA ESTA IGUAL AL COMP ORI
coef1 = 4.0*(alpha*beta)/(1.0 + beta*gamma);
coef2 = 1.0/(1.0 + alpha*gamma) - 2.0*(alpha*beta)/(1.0 + beta*gamma);
dti_b = coef0.*(coef1*aa + coef2*bb);
                   
%% Suma de ambas derivadas
% Smoothing of the topological derivative
pchi = (psi<0); tchi = pdeintrp(p,t,pchi);
% Compute function g: g = -DT (bulk) and g = DT (inclusion)
dtE_m = -tchi.*dte_m + (1-tchi).*dti_m; dtE_m = pdeprtni(p,t,dtE_m);
dtE_b = -tchi.*dte_b + (1-tchi).*dti_b; dtE_b = pdeprtni(p,t,dtE_b);
dtE_shell = dtE_m + dtE_b;
%salida para comp shell
dt=dtE_shell;

end