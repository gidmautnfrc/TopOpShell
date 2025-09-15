%**************************************************************************
% Structural Compliance Topology Optimization with Volume Constraint  
%**************************************************************************
% 
% DESCRIPTION
% Computes the topological derivative and use it together with a 
% level-set domain representation method in the context of structural 
% topology optimization design
%
% HISTORY
% A.A. Novotny     01/2013: code implementation
% A.A. Novotny     02/2020: code updating
%
% MAIN REFERENCE
%
%
%**************************************************************************
clear all; close all; format long e; clc;

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
%     e = mesh.e;
    t = mesh.t;
%     g = mesh.g; 
    area = mesh.A;
    np = mesh.np;
    vol0 = sum(area);
    per_defined= false;
    Per0=1;
    alpha=1;
    params.alpha= alpha;
    perimeter = 1;
    dtper=0;
    
    % topology optimization parameters
    stop = params.stop; 
    kmin = params.kmin; 
    penalty = params.penalty; 
% penalization = params.penalization;
% volfrac = params.volfrac;
% auglag = params.auglag;
% voleps = params.voleps;

%% PLOT 1
% figure(10); clf;
% pdeplot(p,e,t,'xydata',t(4,:),'xystyle','flat','colormap','gray',...
%               'xygrid','off','colorbar','off','title','Geometry'); 
% axis image; axis off;
%% PLOT 1
% figure(10); clf;
% trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),ones(1,length(p)),'LineStyle','none')
% axis off
%% MAS PARAMETROSDE MATRIZ DE MASA UNITARIA PSI Y TCHI

        [~,unitM,~] = assema(p,t,0,1,0); % mass matrix of unity density --> PREGUNTAR AUGUSTO; da distinto al ejemplo del ex1K siendo la misma malla
                                                                           % pero con coord z=0;
        psi = psi/sqrt(dot(unitM*psi,psi)); % level-set function nomalization 
        tchi = pdeintrp(p,t,(psi < 0)); vol = dot(area,tchi); % volume
        M=unitM;
        
%% PLOT         
%         plot initial guess
%         figure(1); clf; set(1,'WindowStyle','docked');
%         pdeplot(p,e,t,'xydata',(psi>0),'xystyle','flat','colormap','gray',...
%                       'xygrid','off','colorbar','off'); axis image; axis off;
%
%% plotear caso normal el U, y sus componentes, F , K y ver las diferencias.

%% PLOT         
% 
        figure('Name','Geometria'); clf; set(1,'WindowStyle','docked');
        trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),'LineStyle','none','FaceColor',[0.5 0.5 0.5]); 
        axis off

%% SOLVE DEL PDE, EN MI CASO A PASAR A SHELLSOLVE
        % hold-all domain
        psi_full = -ones(np,1);              % el vector C en la posicion 10 tiene el auxiliar de gamma para este calculo
%       [U,F] = shellsolve(psi,mesh,aux,pdecoef,bc); % --> no se que calcula inicialmente aca, que calcula con un gamma 1 y un aux, 
                                                     %      y no entiendo como afecta el a y el tgamma al resultado del assema ori
        [U,F] = shellsolve(mesh,pdecoef,matprop,signatures,bc,psi_full);
        comp0 = 0.5*dot(F,U.U_shell);

        % solve linear system
        [U,F] = shellsolve(mesh,pdecoef,matprop,signatures,bc,psi); %% --> supongo que resuelve normal el sistema, con las condiciones y
                                                                     %       todo comun como venia ahciendo
        energy = 0.5*dot(F,U.U_shell); 

        % compute shape function
        sf = energy/comp0 + penalty * vol/vol0;  

%       figure(2); clf; pdesurf(p,t, U(1:np)); 

        k = 1; iter = 0; option = 'null'; 
        gsf = sf; gth = pi; git = iter; gvo = vol;
        
% %% PLOT         
% 
%         figure('Name','otro'); clf; set(1,'WindowStyle','docked');
%         trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),(p(3,:)+U.U_shell(2*np+1:3*np,1)'),'LineStyle','none'); 
%         axis off
% %
        disp = U.U_shell;
        scale = 1;
        disp = disp*scale;
        load_nodes = bc.pNeu(:,1);
        clamped_nodes = bc.pDir(:,1);
% % % %% PLOT SHELL RESULTS
% % % 
%     pplot=p';
%      figure('Name','Resultados','NumberTitle','off')
%      trisurf(t(1:3,:)',pplot(:,1),pplot(:,2.),pplot(:,3),'FaceAlpha',0.7,'FaceColor',[0.58 0.58 0.58]),; xlabel('x'), ylabel('y'), zlabel('z')
% %      legend('Estructura inicial')
%      hold on,
%      trisurf(t(1:3,:)',pplot(:,1)+disp(1:np),pplot(:,2)+disp(np+1:2*np),pplot(:,3)+disp(2*np+1:3*np),'FaceAlpha',0.4,'FaceColor',[0 0 1]);
% %      legend('Estructura deformada')
%      hold on, plot3(pplot(load_nodes,1),pplot(load_nodes,2),pplot(load_nodes,3),'rv')
% %      legend('Carga')
%      hold on,  plot3(pplot(clamped_nodes,1),pplot(clamped_nodes,2),pplot(clamped_nodes,3),'*g','MarkerSize',10)
%      legend('Estructura inicial','Estructura deformada','Carga externa','Empotramiento')
%%
% COLOR MAP DE GRISES
r = [0 0 0];       %# start
w = [.25 .25 .25];    %# middle
b = [0.5 0.5 0.5];       %# end

%# colormap of size 64-by-3, ranging from red -> white -> blue
c1 = zeros(32,3); c2 = zeros(32,3);
for i=1:3
    c1(:,i) = linspace(r(i), w(i), 32);
    c2(:,i) = linspace(w(i), b(i), 32);
end
map = [c1(1:end-1,:);c2];

%%

% % % %% ARRANCA EL while para iterar

while not(strcmp(option,'s'))
        
    iter = iter + 1; 
    if(iter == 200) 
        option = 's';
    end
    
%    dt = topder(U,psi,mesh,matprop,pdecoef); % derivada topologica
    [dt] = tdshell(mesh,U,pdecoef, matprop,signatures,psi);

    if per_defined==true
    [dtper,perimeter] = TopDerPer(mesh,psi,params,per_defined,Per0);
    end

    dt = dt/comp0 + penalty/vol0;
    dt = dt + dtper;
    dt = dt/sqrt(dot(unitM*dt,dt));
    sf = sf + alpha * perimeter/Per0;

    
    cosin = max(min(dot(unitM*dt,psi),1.0),-1.0);
    theta = max(real(acos(cosin)),1.0e-4);

   

    % performs a line-search 
    sfold = sf; psiold = psi; 
    sf = sf + 1.0; k = min(1,1.5*k); % trick
    
 
    while (sf > sfold)
        
        if(k < kmin) 
            theta = stop / 2;

        end                    

        % update level-set function
        psi  = (sin((1-k)*theta)*psiold + sin(k*theta)*dt)./sin(theta);        
        psi = psi/sqrt(dot(unitM*psi,psi));  
        
        % solve linear system
        %[U,F] = pdesolve(psi,mesh,matprop,pdecoef,bc); 
        [U,F] = shellsolve(mesh,pdecoef,matprop,signatures,bc,psi);
        energy = 0.5*dot(F,U.U_shell); 
        if per_defined==true
        [~,perimeter] = TopDerPer(mesh,psi,params,per_defined,Per0);
        end
     
        % update the volume of the bulk phase
        tchi = pdeintrp(p,t,(psi < 0)); vol = dot(area,tchi);     
       
        % compute shape function
        sf = energy/comp0 + penalty * vol/vol0;
        sf = sf + alpha * perimeter/Per0;
        k = k / 2;

    end   
            
    k = k * 2;

    git = [git,iter]; gsf = [gsf,sf]; gth = [gth,theta]; gvo = [gvo,vol];    
    
    disp(['iter   = ', num2str(iter)]);
    disp(['volume = ', num2str(vol),' => ', num2str(vol*100/vol0), '%']);
    disp(['energy = ', num2str(energy)]);
    disp(['shfunc = ', num2str(sf)]); 
    disp(['kappa  = ', num2str(k)]); 
    disp(['theta  = ', num2str(theta*180/pi)]);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
    drawnow
%     figure(3); clf; set(3,'WindowStyle','docked'); %%plot de la topologia
%     pdeplot(p,e,t,'xydata', (psi > 0),'xystyle','flat','colormap','gray',...
%                   'xygrid','off','colorbar','off'); axis image; axis off;           
    
    topplot=1*(psi>0);
    tchilogical=find(tchi>0);
    tplot= t(:,tchilogical);
    tchiplot=-tchi(tchilogical);    
    
    figure(3); clf; set(3,'WindowStyle','docked');
    trisurf(tplot(1:3,:)',p(1,:),p(2,:),p(3,:), tchiplot,'LineStyle','none');
    colormap(map)
    axis equal
    
%     
%     plot = trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:), -tchi,'LineStyle','none');
%     figure(plot); clf; set(3,'WindowStyle','docked');
%     axis equal

    figure(4); clf; set(4,'WindowStyle','docked');
    plot(git,gsf, ':*k'); title('Shape Function');
    
    figure(5); clf; set(5,'WindowStyle','docked');   
    plot(git,gth*180/pi, ':*k'); title('Theta Angle');
    
    figure(6); clf; set(6,'WindowStyle','docked');
    plot(git,gvo/vol0, ':*k'); title('Volume Fraction');   
    
    if or(k < kmin,theta < stop)
        % stop or try a mesh refinement        
        while and(not(strcmp(option,'r')), not(strcmp(option,'s')))
          option = input('\n -> type "r" to remesh or "s" to stop : ', 's');
        end
        
        if (option == 'r')
            cd('examples')
%                 [mesh, pdecoef, matprop, params, bc, psi] = example(mesh,psi,params);
                  [mesh, params,psi, bc, signatures,C] = example(mesh,psi,params);
            cd ..
            % update data            
               % mesh and geometry  parameters
                p = mesh.p; t = mesh.t;
                area = mesh.A; np = mesh.np; vol0 = sum(area);
            
            [~,unitM,~] = assema(p,t,0,1,0);
            psi = psi/sqrt(dot(unitM*psi,psi)); 
            tchi = pdeintrp(p,t,(psi < 0)); 
            vol = dot(area,tchi); % volume
            
%             aux = matprop; aux.gamma = 1.0 ;% hold-all domain
%             [U,F] = pdesolve(psi,mesh,aux,pdecoef,bc); comp0 = dot(F,U);
%             [U,F] = pdesolve(psi,mesh,matprop,pdecoef,bc); energy = dot(F,U); 

             aux = matprop; aux.gamma = 1.0;
            [U,F] = shellsolve(mesh,pdecoef,aux,signatures,bc,psi);
            comp0 = 0.5*dot(F,U.U_shell);
            
            
            [U,F] = shellsolve(mesh,pdecoef,matprop,signatures,bc,psi); 
            energy = 0.5*dot(F,U.U_shell);
            [~,perimeter] = TopDerPer(mesh,psi,params,per_defined,Per0);


            if ~per_defined
                per0= perimeter;
                per_defined=true;
                disp('Se activo la restriccion perimetrica');
            end


            
            sf = energy/comp0 + penalty * vol/vol0;
            sf = sf + alpha * perimeter/Per0;
            k = 1;
            
            option = 'null';
            
        end
    end
end

fig = input(' -> export figures? type "y" or "n"   : ', 's');
if strcmp(fig,'y')
    cd('results')
        print (1, '-dtiff', 'iniguess');
%         print (2, '-dtiff', 'deformed');
        print (3, '-dtiff', 'topology');
        print (4, '-deps',  'shfunc');
        print (5, '-deps',  'theta');
        print (6, '-deps',  'volume');
    cd ..
else
    return;
end