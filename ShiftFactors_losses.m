clc
clearvars
syms b12 b13 b23 g12 g13 g23 d1 d2 d3 Pg1 Pg2 Pg3 Pd2 Pd3 real

Pg = [Pg1; Pg2+Pg3; 0];
Pd = [0; Pd2; Pd3];

Af = [1 0 0;
      1 0 0;
      0 1 0];

At = [0 -1 0;
      0 0 -1
      0 0 -1];

A = Af + At;

y = [g12+1i*b12 g13+1i*b13 g23+1i*b23];

d=[d1;d2;d3];

Bf = diag(imag(y))*A;
Bbus = A'*Bf;

fk = Bf*d;

pl12 = g12*(d1^2 - 2*d1*d2 + d2^2);
pl13 = g13*(d1^2 - 2*d1*d3 + d3^2);
pl23 = g23*(d2^2 - 2*d2*d3 + d3^2);

pl = [pl12;pl13;pl23];

f_f = fk + 0.5*pl;
f_t = fk - 0.5*pl;

Ir = [0 0;
      1 0;
      0 1];

Bbus_r = Ir'*Bbus*Ir;
SF = Bf*Ir*(Bbus_r\Ir');
%%

Bf*Ir/Bbus_r

simplify(SF*Ir)

%%
f12 = b12*(d1-d2);
f13 = b13*(d1-d3);
f23 = b23*(d2-d3);

f12_f = f12 + pl12;
f13_f = f13 + pl13;
f23_f = f23 + pl23;
f12_t = f12 - pl12;
f13_t = f13 - pl13;
f23_t = f23 - pl23;


