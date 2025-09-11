%% TESTEOS
% K
K_uvuv=K_SHELL(1:np*2,1:np*2);
K_tztz=K_SHELL(np*5+1:end,np*5+1:end);
K_tzuv=K_SHELL(1:np*2,np*5+1:end);
K_uvtz=K_SHELL(np*5+1:end,1:np*2);
K_memb=[K_uvuv K_tzuv;
        K_uvtz K_tztz]; 
K_dif=max(abs((K_comp-K_uvuv)));
% Aux

dt_plate=dt; U_plate=U; psi_plate=psi;

U_shell_comp=U.U_shell(2*np+1:5*np,1);
% Errores
[error_rel_psi, nod_psi] = max(abs((psi_plate-psi)./psi_plate));
[error_rel_u, gdl_u] = max(abs((U_plate-U_shell_comp)./U_plate));
[error_rel_u, gdl_u] = max(abs((U_membrane-U_shell_comp)./U_membrane));
[error_rel_dt, nod_dt] = max(abs((dt_plate-dt)./dt_plate));

%% VERIFICACION DE 0

U_plate(gdl_u)==0
U_shell_comp(gdl_u)==0
U_plate(gdl_u)
U_shell_comp(gdl_u)

psi(nod_psi)
psi_plate(nod_psi)

dt_plate(nod_dt)
dt(nod_dt)
%% PLOT nodo        

        figure('Name','TEST'); clf; set(1,'WindowStyle','docked');
        trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),'LineStyle','none','FaceColor',[0.5 0.5 0.5]); 
        hold on, 
        plot3(p(1,nod_psi),p(2,nod_psi),p(3,nod_psi),'*k')
        plot3(p(1,nod_u),p(2,nod_u),p(3,nod_u),'*g')
        plot3(p(1,nod_dt),p(2,nod_dt),p(3,nod_dt),'*r')
        
