function s = getstress_plate(mesh,U, matprop)
%Preliminaries

t = mesh.t;
nt=size(t,2);
Ub=U.Ub;

%Local coordinates

x1 =  mesh.localcoordinates.x1;
x2 =  mesh.localcoordinates.x2;
x3 =  mesh.localcoordinates.x3;
y1 =  mesh.localcoordinates.y1;
y2 =  mesh.localcoordinates.y2;
y3 =  mesh.localcoordinates.y3;
 
x1=x1';
y1=y1';
x2=x2';
y2=y2';
x3=x3';
y3=y3';

%Alpha Matrix
x12 = x1-x2;
y12 = y1-y2;
x23 = x2-x3;
y23 = y2-y3;
x31 = x3-x1;
y31 = y3-y1;

l122 = x12.*x12 + y12.*y12;
l232 = x23.*x23 + y23.*y23;
l312 = x31.*x31 + y31.*y31;

p4 = -6*x23./l232;
p5 = -6*x3 ./l312;
p6 = -6*x12./l122;

t4 = -6*y23./l232;
t5 = -6*y3 ./l312;

q4 = 3*x23.*y23./l232;
q5 = 3*x3.*y3  ./l312;

r4 = 3*(y23.*y23)./l232;
r5 = 3*(y31.*y31)./l312;

TwoArea = x2 .*y3;

  %Material Properties
  
        h0 = matprop.h0; 
        la = matprop.la_dkt; 
        mu = matprop.mu_dkt;     
        nu = la/(2*mu+la); 
        E0 = 2*mu*(1+nu); 

        %Material Properties    
        E1 = E0*h0*h0*h0/(12*(1-nu*nu));
        E2 = nu.*E1;
        E3 = E1;
        E4 = E1.*(1.-nu)./2;


%Elements of stiffness matrix
alfa11 = y3.*p6;
%alfa12 = 0.0;
alfa13 = -4.0.*y3;
alfa14 = -alfa11;
%alfa(1,5) = 0.0;
alfa16 = -2.0.*y3;
%alfa17 = 0.0;
%alfa18 = 0.0;
%alfa19 = 0.0;
alfa21 = alfa14;
%alfa22 = 0.0;
alfa23 = -alfa16;
alfa24 = alfa11;
%alfa25 = 0.0;
alfa26 = -alfa13;
%alfa27 = 0.0;
%alfa28 = 0.0;
%alfa29 = 0.0;
alfa31 = y3.*p5;
alfa32 = -y3.*q5;
alfa33 = y3.*(2.0 - r5);
alfa34 = y3.*p4;
alfa35 = y3.*q4;
alfa36 = y3.*(r4 - 2.0);
alfa37 = -y3.*(p4 + p5);
alfa38 = y3.*(q4 - q5);
alfa39 = y3.*(r4 - r5);
alfa41 = -x2.*t5;
alfa42 = x23 + x2.*r5;
alfa43 = -x2.*q5;
%alfa44 = 0.0;
alfa45 = x3;
%alfa46 = 0.0;
alfa47 = -alfa41;
alfa48 = x2.*(r5 - 1.0);
alfa49 = alfa43;
%alfa51 = 0.0;
alfa52 = x23;
%alfa53 = 0.0;
alfa54 = x2.*t4;
alfa55 = x3 + x2.*r4;
alfa56 = -x2.*q4;
alfa57 = -alfa54;
alfa58 = x2.*(r4 - 1.0);
alfa59 = alfa56;
alfa61 = x23.*t5;
alfa62 = x23.*(1.0 - r5);
alfa63 = x23.*q5;
alfa64 = -x3.*t4;
alfa65 = x3.*(1.0 - r4);
alfa66 = x3.*q4;
alfa67 = -alfa61 -alfa64;
alfa68 = -x23.*r5 - x3.*r4 - x2;
alfa69 = alfa66 + alfa63;
alfa71 = -x3.*p6 - x2.*p5;
alfa72 = -alfa43 + y3;
alfa73 = -4.0.*x23 + x2.*r5;
alfa74 = x3.*p6;
alfa75 = -y3;
alfa76 = 2.0.*x3;
alfa77 = x2.*p5;
alfa78 = -alfa43;
alfa79 = (r5 - 2.0).*x2;
alfa81 = -x23.*p6;
alfa82 = y3;
alfa83 = 2.0.*x23;
alfa84 = -alfa81 + x2.*p4;
alfa85 = -y3 - alfa56;
alfa86 = -4.0.*x3 + x2.*r4;
alfa87 = -x2.*p4;
alfa88 = -alfa56;
alfa89 = (r4 - 2.0).*x2;
alfa91 = x23.*p5 + y3.*t5;
alfa92 = -alfa63 + (1.0 - r5).*y3;
alfa93 = (2.0 - r5).*x23 - alfa32;
alfa94 = -x3.*p4 + y3.*t4;
alfa95 = (r4 - 1.0).*y3 - alfa66;
alfa96 = (2.0 - r4).*x3 - alfa35;
alfa97 = -x23.*p5 + x3.*p4 - (t4 + t5).*y3;
alfa98 = -alfa63 - alfa66 + (r4 - r5).*y3;
alfa99 = -x23.*r5 - x3.*r4 + 4.0*x2 + (q5 - q4).*y3;

% Vetor deslocamento

U1 = Ub(1,:)';       %w1
U2 = Ub(4,:)';       %titax1
U3 = Ub(7,:)';       %titay1
U4 = Ub(2,:)';       %w2
U5 = Ub(5,:)';       %titax2
U6 = Ub(8,:)';       %titay2
U7 = Ub(3,:)';       %w3
U8 = Ub(6,:)';       %titax3
U9 = Ub(9,:)';       %titay3

Cos=1;
Sen=0;

StrL(1,:) = (1/3./TwoArea).*((((alfa11 + alfa21 + alfa31).*U1 + (alfa14 + alfa24 + alfa34).*U4 + alfa37.*U7) ...
    + (alfa32.*U2 + (alfa13 + alfa23 + alfa33).*U3 + alfa35.*U5 + (alfa16 + alfa26 + alfa36).*U6 + alfa38.*U8 + alfa39.*U9).*Cos - ((alfa13 + alfa23 + alfa33).*U2 - alfa32.*U3 + (alfa16 + alfa26 + alfa36).*U5 - alfa35.*U6 + alfa39.*U8 - alfa38.*U9).*Sen).*E1 ...
    + (((alfa41 + alfa61).*U1 + (alfa54 + alfa64).*U4 + (alfa47 + alfa57 + alfa67).*U7) ...
    + ((alfa42 + alfa52 + alfa62).*U2 + (alfa43 + alfa63).*U3 + (alfa45 + alfa55 + alfa65).*U5 + (alfa56 + alfa66).*U6 + (alfa48 + alfa58 + alfa68).*U8 + (alfa49 + alfa59 + alfa69).*U9).*Cos - ((alfa43 + alfa63).*U2 - (alfa42 + alfa52 + alfa62).*U3 + (alfa56 + alfa66).*U5 - (alfa45 + alfa55 + alfa65).*U6 + (alfa49 + alfa59 + alfa69).*U8 - (alfa48 + alfa58 + alfa68).*U9).*Sen).*E2);

StrL(2,:) = (1/3./TwoArea).*((((alfa11 + alfa21 + alfa31).*U1 + (alfa14 + alfa24 + alfa34).*U4 + alfa37.*U7) ...
    + (alfa32.*U2 + (alfa13 + alfa23 + alfa33).*U3 + alfa35.*U5 + (alfa16 + alfa26 + alfa36).*U6 + alfa38.*U8 + alfa39.*U9).*Cos - ((alfa13 + alfa23 + alfa33).*U2 - alfa32.*U3 + (alfa16 + alfa26 +alfa36).*U5 - alfa35.*U6 + alfa39.*U8 - alfa38.*U9).*Sen).*E2 ...
    + (((alfa41 + alfa61).*U1 + (alfa54 + alfa64).*U4 + (alfa47 + alfa57 + alfa67).*U7) ...
    + ((alfa42 + alfa52 + alfa62).*U2 + (alfa43 + alfa63).*U3 + (alfa45 + alfa55 + alfa65).*U5 + (alfa56 + alfa66).*U6 + (alfa48 + alfa58 + alfa68).*U8 + (alfa49 + alfa59 + alfa69).*U9).*Cos - ((alfa43 + alfa63).*U2 - (alfa42 + alfa52 + alfa62).*U3 + (alfa56 + alfa66).*U5 - (alfa45 + alfa55 + alfa65).*U6 + (alfa49 + alfa59 + alfa69).*U8 - (alfa48 + alfa58 + alfa68).*U9).*Sen).*E3);

StrL(3,:) = (1/3./TwoArea).*((((alfa71 + alfa81 + alfa91).*U1 + (alfa74 + alfa84 + alfa94).*U4 + (alfa77 + alfa87 + alfa97).*U7) ...
    + ((alfa72 + alfa82 + alfa92).*U2 + (alfa73 + alfa83 + alfa93).*U3 + (alfa75 + alfa85 + alfa95).*U5 + (alfa76 + alfa86 + alfa96).*U6 + (alfa78 + alfa88 + alfa98).*U8 + (alfa79 + alfa89 + alfa99).*U9).*Cos - ((alfa73 + alfa83 + alfa93).*U2 - (alfa72 + alfa82 + alfa92).*U3 + (alfa76 + alfa86 + alfa96).*U5 - (alfa75 + alfa85 + alfa95).*U6 + (alfa79 + alfa89 + alfa99).*U8 - (alfa78 + alfa88 + alfa98).*U9).*Sen).*E4);

% Rotaciona Deformacao e tensao
cc = Cos.*Cos;
ss = Sen.*Sen;
sc = Sen.*Cos;

% Global s = [sxx, sxy, syy]
s=zeros(3,nt);
s(1,:) = (cc.*StrL(1,:) + ss.*StrL(2,:) - 2.0.*sc.*StrL(3,:));
s(2,:) = (sc.*StrL(1,:) - sc.*StrL(2,:) + (cc - ss).*StrL(3,:));
s(3,:) = (ss.*StrL(1,:) + cc.*StrL(2,:) + 2.0.*sc.*StrL(3,:));

end