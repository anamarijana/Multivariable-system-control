clear; 
close all;
clc

%% linearizovani sistem
kapa_G=nelin_to_lin;
kapa_G_tf = tf(kapa_G);
%% Nominalne vrednosti

kw=1e-14;      %mol^2/(dm^3)^2 
Ca = 1e-6;     %mol/(dm^3)
Cb = 1e-6;     %mol/(dm^3)
V = 30;        %dm^3
kv = 0.01;     %l/s
Fae = 1/60;    %l/s 
Fbe = 2/60;    %l/s

xae = Fae*Ca/(Fae+Fbe);
xbe = Fbe*Cb/(Fae+Fbe+kv); 

H_plus_e = sqrt(power((xbe-xae),2)/4+kw) - (xbe-xae)/2; 
pHe = -log10(H_plus_e);

%% Statičko skaliranje objekta upravljanja
pHmax =  0.0596; %poteklo iz odziva linearni_nelinearni_matlab.m za data maksimalna upravljanja
xamax = 2.0819e-08; %takodje
Famax = 0.1*Fae;
Fbmax = 0.1*Fbe;
Dy = [pHmax 0; 0 xamax];
Du = [Famax 0; 0 Fbmax];

G = inv(Dy)*kapa_G*Du;
G_tf = tf(G);
figure;
step(G)

delta_pH = pHmax/pHe;
delta_xa = xamax/xae;
%% Zanemarivanje
C = [0 1; 1 0];
Gstar = G*inv(C);
Gstar_tf = zpk(tf(Gstar));
Gstar11 = Gstar_tf(1,1);
Gstar22 = Gstar_tf(2,2);
Gstar21 = Gstar_tf(2,1);
Gstar12 = Gstar_tf(1,2);

figure;
step(Gstar)

%% 

w1_design = 0.0017;
Fpf_design=90*pi/180;
s = tf('s');

Kstar11 = minreal(w1_design/s*inv(Gstar11));
Kstar22 = minreal(w1_design/s*inv(Gstar22));

%wpk = 10*w1_design;
%Kstar11 = Kstar11/(s/wpk+1);
%Kstar22 = Kstar22/(s/wpk+1);

% Ti1 = 1/w1_design*tan(Fpf_design-pi/2-unwrap(angle(freqresp(Gstar11,w1_design))));
% Kc1 = 1/abs(freqresp((1+1/Ti1/s)*Gstar11,w1_design));
% Kstar11 = Kc1*(1 + 1/Ti1/s);
% figure; margin(Kstar11*Gstar11);
% 
% Ti2 = 1/w1_design*tan(Fpf_design-pi/2-unwrap(angle(freqresp(Gstar22,w1_design))));
% Kc2 = 1/abs(freqresp((1+1/Ti2/s)*Gstar22,w1_design));
% Kstar22 = Kc2*(1 + 1/Ti2/s);
% figure; margin(Kstar22*Gstar22);
%%
Wstar11 = Gstar11*Kstar11;
Wstar12= Gstar12*Kstar11;
Wstar22 = Gstar22*Kstar22;
Wstar21 = Gstar21*Kstar22;
%%
Kstar = [Kstar11 0; 0 Kstar22];
K   = inv(C)*Kstar;

W = minreal(G_tf*K);
W_tf = zpk(minreal(tf(W)));
Wpoles = pole(W);
detW = minreal(W_tf(1,1)*W_tf(2,2)-W_tf(2,1)*W_tf(1,2));

%%
T = minreal(G_tf*K*inv(eye(2)+G_tf*K));

figure;
step(W_tf);
title('NESTABILNO otvorena sprega dijagonalizovan objekat')

%%
[response, timearray] = step(T);

response11 = response(:,1,1); 
response12 = response(:,2,1); 
response21 = response(:,1,2);
response22 = response(:,2,2);

figure;
subplot(2,2,1)
plot(timearray,response11,'k','LineWidth',2);
ylabel('pH [au]');xlabel('vreme [s]')
xlim([min(timearray),max(timearray)])
grid
title('\rm Fb')
subplot(2,2,2)
plot(timearray,response21,'k','LineWidth',2);
ylabel('pH [au]');xlabel('vreme [s]')
xlim([min(timearray),max(timearray)])
title('\rm Fa')
grid
subplot(2,2,3)
plot(timearray,response12,'k','LineWidth',2);
ylabel('xa [mol/l]');xlabel('vreme [s]')

xlim([min(timearray),max(timearray)])
grid
title('\rm Fb')
subplot(2,2,4)
plot(timearray,response22,'k','LineWidth',2);
ylabel('xa [mol/l]');xlabel('vreme [s]')
xlim([min(timearray),max(timearray)])
grid
title('\rm Fa')
sgtitle({'Odskočni odziv zatvorene sprege',' decentralizovan kontroler',' zanemarivanje sprežnih elemenata'})%,'fontweight','bold','fontsize',14);

%% Singularne karakteristike otvorena sprega, informativno
w=logspace(-5,-1,500);

W=minreal(G*K);           %otvorena sprega
sv_W = sigma(W,w);          %singularne karakteristike otvorene sprege
sv_G = sigma(G,w);

w=logspace(-5,-1,500);
figure;
subplot(2,1,1);
semilogx(w,20*log10(sv_G(1,:)),'b','LineWidth',2);
hold on;
semilogx(w,20*log10(sv_G(2,:)),'r','LineWidth',2);
hold off;
set(gca, 'XScale', 'log');
xlabel('\omega [rad/s]'); ylabel('\sigma(G(j\omega))_d_B'); 

subplot(2,1,2); 
semilogx(w,20*log10(sv_W(1,:)),'b','LineWidth',2);
hold on;
semilogx(w,20*log10(sv_W(2,:)),'r','LineWidth',2);
hold off;
set(gca, 'XScale', 'log');
xlabel('\omega [rad/s]'); ylabel('\sigma(W(j\omega))_d_B'); 
sgtitle('Singularne karakteristike otvorene sprege, bez i uz kontroler')

%% Signularne karakteristike zatvorene sprege, uz kontroler

T = W*(eye(2)+W)^(-1); % komplementarna osetljivost (spregnuti prenos)
sv_T = sigma(T,w); %singularne karakteristike kompl osetljivost

S = eye(2) - T; %osetljivost (spregnuti prenos)
sv_S=sigma(S,w); %singularne karakteristike osetljivost

figure;
subplot(2,1,1);
semilogx(w,20*log10(sv_S(1,:)),'b','LineWidth',2);
hold on;
semilogx(w,20*log10(sv_S(2,:)),'r','LineWidth',2);
hold off;
set(gca, 'XScale', 'log');
xlabel('\omega [rad/s]'); ylabel('\sigma(S(j\omega))_d_B'); 

subplot(2,1,2); 
semilogx(w,20*log10(sv_T(1,:)),'b','LineWidth',2);
hold on;
semilogx(w,20*log10(sv_T(2,:)),'r','LineWidth',2);
hold off;
set(gca, 'XScale', 'log');
xlabel('\omega [rad/s]'); ylabel('\sigma(T(j\omega))_d_B'); 
sgtitle({'Singularne karakteristike zatvorene sprege',' decentralizovan kontroler',' zanemarivanje sprežnih elemenata'})

%% Frekvencijski pokazatelji performansi

gornja_sv_S = sv_S(1,:);
gornja_sv_T = sv_T(1,:);

w0 = w(find(gornja_sv_S >= 1/sqrt(2),1))
MS = max(gornja_sv_S)
MT = max(gornja_sv_T)





