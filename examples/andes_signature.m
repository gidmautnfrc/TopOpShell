function [signature]=andes_signature(firm,nu)
%% Save the signature of andes element
% Entry('',C) = 'OPT' -> Optimal Andes Signature
%             '3I' -> Optimal Andes Signature  
%             '3M' -> Optimal Andes Signature
%             'LS' -> Optimal Andes Signature
% 
% Out = Struct with the signatures
%% ANDES SIGNATURES

    switch(firm)
        case 'OPT'
                %% Firma ANDES Optimo
                  signature.Alphab =  1.5;
                  signature.Beta0 =  0.5*(1-4*nu^2);
                  signature.Beta1 =  1;
                  signature.Beta2 =  2;
                  signature.Beta3 =  1;
                  signature.Beta4 =  0;
                  signature.Beta5 =  1;
                  signature.Beta6 =  -1;
                  signature.Beta7 =  -1;
                  signature.Beta8 =  -1;
                  signature.Beta9 =  -2;
        case '3I'
                %% Firma ALLMAN 3I
                  signature.Alphab =  1;
                  signature.Beta0 =  4/9;
                  signature.Beta1 =  1/12;
                  signature.Beta2 =  5/12;
                  signature.Beta3 =  1/2;
                  signature.Beta4 =  0;
                  signature.Beta5 =  1/3;
                  signature.Beta6 =  -1/3;
                  signature.Beta7 =  -1/12;
                  signature.Beta8 =  -1/2;
                  signature.Beta9 =  -5/12;
        case '3M'
                %% Firma ALLMAN 3M
                  signature.Alphab =  1;
                  signature.Beta0 =  4/9;
                  signature.Beta1 =  1/4;
                  signature.Beta2 =  5/4;
                  signature.Beta3 =  3/2;
                  signature.Beta4 =  0;
                  signature.Beta5 =  1;
                  signature.Beta6 =  -1;
                  signature.Beta7 =  -1/4;
                  signature.Beta8 =  -3/2;
                  signature.Beta9 =  -5/4;
        case 'LS'
                %% Firma ALLMAN LS
                  signature.Alphab =  1;
                  signature.Beta0 =  4/9;
                  signature.Beta1 =  3/20;
                  signature.Beta2 =  3/4;
                  signature.Beta3 =  9/10;
                  signature.Beta4 =  0;
                  signature.Beta5 =  3/5;
                  signature.Beta6 =  -3/5;
                  signature.Beta7 =  -3/20;
                  signature.Beta8 =  -9/10;
                  signature.Beta9 =  -3/4;
        otherwise  
                 disp('Do not match with any signature')
    end
    %% Coeficiente para las deformaciones
    % Equivalente a CST
        signature.dseta1=1/3;
        signature.dseta2=1/3;
        signature.dseta3=1/3;
    
end