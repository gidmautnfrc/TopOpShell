function  [U] = disp_filter(mesh, U_shell);
% variable entrada, mesh y desplazamiento de cascaras como:
%
%     U = [u1
%          .
%          .   
%          un
%          v1
%          . 
%          .          
%          vn
%          w1
%          . 
%          .          
%          wn
%        thetax1
%          . 
%          .          
%        thetaxn
%        thetay1
%          . 
%          .          
%        thetan%
%        thetaz1
%          . 
%          .          
%        thetazn];
%
% y la salida es el U filtrado a los desplazamientos de bending (Ub) y los
% desplazameintos de membrane (Um) en coordenadas locales a cada elemento
% finito
%
% U.Um = [ u1
%          .
%          .   
%          un
%          v1
%          . 
%          .          
%          vn
%        thetaz
%          . 
%          .          
%        thetazn];
%
%  U.Ub =  [ w1
%            .
%            .   
%            wn
%         thetax1
%            . 
%            .          
%         thetaxn
%         thetay
%            . 
%            .          
%         thetayn];
%
%%    
 U.U_shell=U_shell;
%%
    
    lambdaxgxl = mesh.lambda.lambdaxgxl; 
    lambdaxgyl = mesh.lambda.lambdaxgyl; 
    lambdaxgzl = mesh.lambda.lambdaxgzl;

    lambdaygxl = mesh.lambda.lambdaygxl;
    lambdaygyl = mesh.lambda.lambdaygyl;
    lambdaygzl = mesh.lambda.lambdaygzl;

    lambdazgxl = mesh.lambda.lambdazgxl;
    lambdazgyl = mesh.lambda.lambdazgyl;
    lambdazgzl = mesh.lambda.lambdazgzl;
    
    t=mesh.t;
    np=mesh.np;
%%

%     Ug1=[U_shell(t(1,:),1)';U_shell(t(1,:)+np,1)';U_shell(t(1,:)+2*np,1)'];
%     Ug2=[U_shell(t(2,:),1)';U_shell(t(2,:)+np,1)';U_shell(t(2,:)+2*np,1)'];
%     Ug3=[U_shell(t(3,:),1)';U_shell(t(3,:)+np,1)';U_shell(t(3,:)+2*np,1)'];
% 
%     Thetag1=[U_shell(t(1,:)+3*np,1)';U_shell(t(1,:)+4*np,1)';U_shell(t(1,:)+5*np,1)'];
%     Thetag2=[U_shell(t(2,:)+3*np,1)';U_shell(t(2,:)+4*np,1)';U_shell(t(2,:)+5*np,1)'];
%     Thetag3=[U_shell(t(3,:)+3*np,1)';U_shell(t(3,:)+4*np,1)';U_shell(t(3,:)+5*np,1)'];

%% Separacion de desplazamientos de U globales a U globales por nodo de elemento por cada elemento
%     Ug_Shell=[Ug1;Thetag1;Ug2;Thetag2;Ug3;Thetag3];
    
    Ug1=U_shell(t(1,:))';
    Vg1=U_shell(t(1,:)+np)';
    Wg1=U_shell(t(1,:)+2*np)';
    Thetax1=U_shell(t(1,:)+3*np)';
    Thetay1=U_shell(t(1,:)+4*np)';
    Thetaz1=U_shell(t(1,:)+5*np)';
    Ug2=U_shell(t(2,:))';
    Vg2=U_shell(t(2,:)+np)';
    Wg2=U_shell(t(2,:)+2*np)';
    Thetax2=U_shell(t(2,:)+3*np)';
    Thetay2=U_shell(t(2,:)+4*np)';
    Thetaz2=U_shell(t(2,:)+5*np)';
    Ug3=U_shell(t(3,:))';
    Vg3=U_shell(t(3,:)+np)';
    Wg3=U_shell(t(3,:)+2*np)';
    Thetax3=U_shell(t(3,:)+3*np)';
    Thetay3=U_shell(t(3,:)+4*np)';
    Thetaz3=U_shell(t(3,:)+5*np)';
    
%% Rotacion global a local

    [Ul1]=lambdaxgxl.*Ug1+lambdaxgyl.*Vg1+lambdaxgzl.*Wg1;
    [Ul2]=lambdaygxl.*Ug1+lambdaygyl.*Vg1+lambdaygzl.*Wg1;
    [Ul3]=lambdazgxl.*Ug1+lambdazgyl.*Vg1+lambdazgzl.*Wg1;
    [Ul4]=lambdaxgxl.*Thetax1+lambdaxgyl.*Thetay1+lambdaxgzl.*Thetaz1;
    [Ul5]=lambdaygxl.*Thetax1+lambdaygyl.*Thetay1+lambdaygzl.*Thetaz1;
    [Ul6]=lambdazgxl.*Thetax1+lambdazgyl.*Thetay1+lambdazgzl.*Thetaz1;
    [Ul7]=lambdaxgxl.*Ug2+lambdaxgyl.*Vg2+lambdaxgzl.*Wg2;
    [Ul8]=lambdaygxl.*Ug2+lambdaygyl.*Vg2+lambdaygzl.*Wg2;
    [Ul9]=lambdazgxl.*Ug2+lambdazgyl.*Vg2+lambdazgzl.*Wg2;
    [Ul10]=lambdaxgxl.*Thetax2+lambdaxgyl.*Thetay2+lambdaxgzl.*Thetaz2;
    [Ul11]=lambdaygxl.*Thetax2+lambdaygyl.*Thetay2+lambdaygzl.*Thetaz2;
    [Ul12]=lambdazgxl.*Thetax2+lambdazgyl.*Thetay2+lambdazgzl.*Thetaz2;
    [Ul13]=lambdaxgxl.*Ug3+lambdaxgyl.*Vg3+lambdaxgzl.*Wg3;
    [Ul14]=lambdaygxl.*Ug3+lambdaygyl.*Vg3+lambdaygzl.*Wg3;
    [Ul15]=lambdazgxl.*Ug3+lambdazgyl.*Vg3+lambdazgzl.*Wg3;
    [Ul16]=lambdaxgxl.*Thetax3+lambdaxgyl.*Thetay3+lambdaxgzl.*Thetaz3;
    [Ul17]=lambdaygxl.*Thetax3+lambdaygyl.*Thetay3+lambdaygzl.*Thetaz3;
    [Ul18]=lambdazgxl.*Thetax3+lambdazgyl.*Thetay3+lambdazgzl.*Thetaz3;
%%

Um=[Ul1; Ul7; Ul13; Ul2; Ul8; Ul14; Ul6; Ul12; Ul18];
Ub=[Ul3; Ul9; Ul15; Ul4; Ul10; Ul16; Ul5; Ul11; Ul17];

    U.Um=Um;
    U.Ub=Ub;
    
end