clear all
close all
clc


%% Parametri

tfin=5e4;     % sec

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


%% Ponašanje oko nominalnog režima

xa0 = xae; 
xb0 = xbe;

pHstep =  0.0596; %poteklo iz odziva linearni_nelinearni_matlab.m za data maksimalna upravljanja
xastep = 2.0819e-08; %takodje

tu1_step1 = 1000;
tu1_step2 = 3000;

tu2_step1 = 5000;
tu2_step2 = 7000;

%% Decentralizovan kontroler: zanemarivanje sprežnih elemenata
s = tf('s');

skaliran_Kstar_dec_z = zpk([0.80709*(s+0.002)*(s+0.001667)/s/(s+0.00181) 0; 0 0.95559*(s+0.001667)/s]);

pHmax =  0.0596; 
xamax = 2.0819e-08;
Famax = 0.1*Fae;
Fbmax = 0.1*Fbe;

Dy = [pHmax 0; 0 xamax];
Du = [Famax 0; 0 Fbmax];

C = [0 1; 1 0];

Kstar_dec_z = Du*skaliran_Kstar_dec_z*inv(Dy);
K_dec_z   = inv(C)*Kstar_dec_z;
%% Decentralizovan kontroler: statičko dekuplovanje

%% Decentralizovan kontroler: dekuplovanje na učestanosti propusnog opsega


%% Pokretanje simulacije
sim('kontrola_nelinearni_sim.slx');

%% nelinearni sistem

min_xa_out_plot = 0;% min(xa_out);
max_xa_out_plot = max(xa_out);
delta_xa_out = max_xa_out_plot-xa_out(1);

min_pH_out_plot = 0; %min(pH_out);
max_pH_out_plot = max(pH_out);
delta_pH_out = max_pH_out_plot-pH_out(1);

min_Fa_out_plot = 0; %min(Fa_out);
max_Fa_out_plot = max(Fa_out);

min_Fb_out_plot = 0; % min(Fb_out);
max_Fb_out_plot = max(Fb_out);

figure;
hold all;
%ylim([min_pH_out_plot,max_pH_out_plot]);
plot(t_out,pH_ref_out,'k--','LineWidth',2)
plot(t_out,pH_out,'k','LineWidth',2)
grid
xlabel('vreme [s]')
ylabel('pH [au]')
title('pH vrednost tečnosti u tenku')
legend('referentna vrednost','decentralizovan, zanemarivanje')
%%
figure;
hold all;
%ylim([min_xa_out_plot,max_xa_out_plot])
plot(t_out,xa_ref_out,'k--','LineWidth',2)
plot(t_out,xa_out,'k','LineWidth',2)
grid
xlabel('vreme [s]')
ylabel('xa [mol/l]')
title('Koncentracija tečnosti u tenku')
legend('referentna vrednost','decentralizovan, zanemarivanje')
%%
figure;
ylim([min_Fa_out_plot,max_Fa_out_plot])
plot(t_out,Fa_out,'k','LineWidth',2)
grid
ylabel('Fa [l/s]')
xlabel('vreme [s]')
title('Protok kiseline u reaktor')
legend('decentralizovan, zanemarivanje')

figure;
ylim([min_Fb_out_plot,max_Fb_out_plot]);
plot(t_out,Fb_out,'k','LineWidth',2)
grid
ylabel('Fb [l/s]')
xlabel('vreme [s]')
title('Protok baze u reaktor')
legend('decentralizovan, zanemarivanje')

