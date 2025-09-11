P = p + 0.001;
[~,index] = min( ((P(1,1)-p(1,:)).^2 + (P(2,1)-p(2,:)).^2).^0.5 );
dummy = ((P(1,:)-p(1,:)).^2 + (P(2,:)-p(2,:)).^2).^0.5 ;
dummy = ((P(1,:)'-p(1,:)).^2 + (P(2,:)'-p(2,:)).^2).^0.5 ;
Requested 65921x65921 (32.4GB) array exceeds maximum array size preference. Creation of arrays greater than this limit may take a long time and cause MATLAB to
become unresponsive. See array size limit or preference panel for more information.
 
dummy = ((P(1,1:10)'-p(1,1:5)).^2 + (P(2,1:10)'-p(2,1:5)).^2).^0.5 ;
[~,index] = min( dummy ,1);
Error using min
MIN with two matrices to compare and two output arguments is not supported.
 
[~,index] = min( dummy ,[],1);
[~,index] = min( dummy ,[],2);
dummy = (psi<0);
charfun = pdeintrp(p,t,dummy);
unique(charfun)

ans =

                         0     3.333333333333333e-01     6.666666666666666e-01     1.000000000000000e+00

