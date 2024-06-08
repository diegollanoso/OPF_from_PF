clc
clearvars

g12=1.9231;
b12=9.6154;
g13=1.1765;
b13=4.7059;
g23=0.9174;
b23=3.0581;


Af = [1 0 0;
      1 0 0;
      0 1 0];

At = [0 -1 0;
      0 0 -1
      0 0 -1];

A = Af + At;

y = [g12+1i*b12 g13+1i*b13 g23+1i*b23];

Bf = diag(imag(y))*A;
Bbus = A'*Bf;


Ir = [0 0;
      1 0;
      0 1];

Bbus_r = Ir'*Bbus*Ir;
SF = Bf*Ir*(Bbus_r\Ir');


%%

SF*Ir

Bf*Ir/(Bbus_r)