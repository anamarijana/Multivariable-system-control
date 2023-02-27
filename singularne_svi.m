clear all;
close all;
clc;

w=logspace(-5,-1,500);

%% linearizovani sistem
kapa_G=nelin_to_lin;
kapa_G_tf = tf(kapa_G);


%% Statičko skaliranje objekta upravljanja
Fae = 1/60;    %l/s 
Fbe = 2/60;    %l/s
pHmax =  0.0596; %poteklo iz odziva linearni_nelinearni_matlab.m za data maksimalna upravljanja
xamax = 2.0819e-08; %takodje
Famax = 0.1*Fae;
Fbmax = 0.1*Fbe;
Dy = [pHmax 0; 0 xamax];
Du = [Famax 0; 0 Fbmax];

G = inv(Dy)*kapa_G*Du;
G_tf = tf(G);

%% Zanemarivanje

C = [0 1; 1 0];
Gstar = G*inv(C);
Gstar_tf = zpk(tf(Gstar));
Gstar11 = Gstar_tf(1,1);
Gstar22 = Gstar_tf(2,2);
Gstar21 = Gstar_tf(2,1);
Gstar12 = Gstar_tf(1,2);

w1_design = 0.0017;
Fpf_design=90*pi/180;
s = tf('s');

Kstar11 = minreal(w1_design/s*inv(Gstar11));
Kstar22 = minreal(w1_design/s*inv(Gstar22));

Kstar = [Kstar11 0; 0 Kstar22];
K_dec_z = inv(C)*Kstar;

%% Staticki

W1_d0 =  inv(dcgain(G));
Gstar_d0 = G*W1_d0;     % Sistem dekuplovan u ustaljenom stanju, kako ga vidi kontroler

Gstar_tf = zpk(tf(Gstar_d0));
Gstar11 = Gstar_tf(1,1);
Gstar22 = Gstar_tf(2,2);

w1_design = 0.0017;
Fpf_design=100*pi/180;
s = tf('s');

%zadajem kakav ocu prenos 
%invertor dinamike
Kstar11 = minreal(w1_design/s*inv(Gstar11));
Kstar22 = minreal(w1_design/s*inv(Gstar22));

Kstar = [Kstar11 0; 0 Kstar22];
K_dec_stat = W1_d0*Kstar;

%% Na ucestanosti
W1_dw =  real_inverse(freqresp(G,w1_design));
Gstar_dw = G*W1_dw;     % Sistem dekuplovan na ucestanosti, kako ga vidi kontroler

Gstar_tf = zpk(tf(Gstar_dw));
Gstar11 = Gstar_tf(1,1);
Gstar22 = Gstar_tf(2,2);

Kstar11 = minreal(w1_design/s*inv(Gstar11));
Kstar22 = minreal(w1_design/s*inv(Gstar22));

Kstar = [Kstar11 0; 0 Kstar22];
K_dec_w0 = W1_dw*Kstar;

%% Hinf

Mswc = 1.2; w0wc = 0.0017; S0wc = 10^(-5);
ws = (s/Mswc+w0wc)/(s+w0wc*S0wc);
Ws1 = [ws 0;0 ws];
Wks = ss((1/16.875)*eye(2));

[K_hinf, N, Hinf] = mixsyn(G,Ws1,Wks,[],'DISPLAY','on');

%% potpuna inverzija

w0 = 0.0017; % zeljeni propusni opseg

W = inv(G); % dekupler
K = eye(2)*(w0/s); % decentralizovani kontroler po kanalima

K_inv = zpk(W*K);
%% Sve singularne karakteristike
S_dec_z = sigma(inv(eye(2) + G*K_dec_z),w);
T_dec_z = sigma(G*K_dec_z*inv(eye(2) + G*K_dec_z),w);

S_dec_stat = sigma(inv(eye(2) + G*K_dec_stat),w);
T_dec_stat = sigma(G*K_dec_stat*inv(eye(2) + G*K_dec_stat),w);

S_dec_w0 = sigma(inv(eye(2) + G*K_dec_w0),w);
T_dec_w0 = sigma(G*K_dec_w0*inv(eye(2) + G*K_dec_w0),w);

S_hinf = sigma(inv(eye(2) + G*K_hinf),w);
T_hinf = sigma(G*K_hinf*inv(eye(2) + G*K_hinf),w);

S_inv = sigma(inv(eye(2) + G*K_inv),w);
T_inv = sigma(G*K_inv*inv(eye(2) + G*K_inv),w);

figure;
subplot(2,1,1);
hold on;
semilogx(w,20*log10(S_dec_z(1,:)),'r');
semilogx(w,20*log10(S_dec_stat(1,:)),'g');
semilogx(w,20*log10(S_dec_w0(1,:)),'b');
semilogx(w,20*log10(S_hinf(1,:)),'c');
semilogx(w,20*log10(S_inv(1,:)),'m');
semilogx(w,20*log10(S_dec_z(2,:)),'r');
semilogx(w,20*log10(S_dec_stat(2,:)),'g');
semilogx(w,20*log10(S_dec_w0(2,:)),'b');
semilogx(w,20*log10(S_hinf(2,:)),'c');
semilogx(w,20*log10(S_inv(2,:)),'m');
hold off;
set(gca, 'XScale', 'log');
xlabel('\omega [rad/s]'); ylabel('\sigma(S(j\omega))_d_B');
legend('dec, zan','dec stat','dec w0','hinf','inv')

subplot(2,1,2);
hold on;
semilogx(w,20*log10(T_dec_z(1,:)),'r');
semilogx(w,20*log10(T_dec_z(2,:)),'r');
semilogx(w,20*log10(T_dec_stat(1,:)),'g');
semilogx(w,20*log10(T_dec_stat(2,:)),'g');
semilogx(w,20*log10(T_dec_w0(1,:)),'b');
semilogx(w,20*log10(T_dec_w0(2,:)),'b');
semilogx(w,20*log10(T_hinf(1,:)),'c');
semilogx(w,20*log10(T_hinf(2,:)),'c');
semilogx(w,20*log10(T_inv(1,:)),'m');
semilogx(w,20*log10(T_inv(2,:)),'m');
hold off;
set(gca, 'XScale', 'log');
xlabel('\omega [rad/s]'); ylabel('\sigma(T(j\omega))_d_B'); 
sgtitle({'Singularne karakteristike zatvorene sprege'})