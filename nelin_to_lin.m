function linsys=nelin_to_lin()
kw=1e-14;      %mol^2/(dm^3)^2 
Ca = 1e-6;     %mol/(dm^3)
Cb = 1e-6;     %mol/(dm^3)
V = 30;        %dm^3
kv = 0.01;     %l/s

%% Nominalni režim - biram nominalno upravljanje Fa i Fb 

Fae = 1/60;    %l/s 
Fbe = 2/60;    %l/s

xae = Fae*Ca/(Fae+Fbe);
xbe = Fbe*Cb/(Fae+Fbe+kv); 

H_plus_e = sqrt(power((xbe-xae),2)/4+kw) - (xbe-xae)/2; % od sada smatram da mi je ovo izlaz
pHe = -log10(H_plus_e);
%% Linearizacija oko nominalnog režima

syms s_xa s_xb s_Fa s_Fb s_eta

%jednačine stanja
f1 = (s_Fa*Ca - (s_Fa + s_Fb)*s_xa)/V;
f2 = (s_Fb*Cb - (s_Fa + s_Fb + kv)*s_xb)/V;

%jednačina izlaza
g = sqrt(power((s_xb-s_xa),2)/4+kw) - (s_xb-s_xa)/2;
h = -log10(g);


A = [diff(f1,s_xa) diff(f1,s_xb);
    diff(f2,s_xa) diff(f2,s_xb)];

B = [diff(f1,s_Fa) diff(f1,s_Fb);
    diff(f2,s_Fa) diff(f2,s_Fb)];

C1 =[diff(g,s_xa) diff(g,s_xb); %gornji: koncentracija vodonika donji: xa 
    1 0];
C2 =[diff(h,s_xa) diff(h,s_xb); %gornji: pH donji: xa 
    1 0];

D = [0 0;0 0];

s_xa = xae; s_xb = xbe; 
s_Fa = Fae; s_Fb = Fbe;
A = eval(A); B = eval(B); 
C = eval(C2);
% C = eval(C2);
%% linearni sistem u prostoru stanja
linsys = minreal(ss(A, B, C, D));
