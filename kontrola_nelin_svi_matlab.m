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

%pHstep =  0.0596; %poteklo iz odziva linearni_nelinearni_matlab.m za data maksimalna upravljanja
pHstep =  0.001*0.0596; %za direkcionalnost minimalnog pojacanja
%xastep = 2.0819e-08; %takodje
xastep = pHstep*0.0003816; %za direkcionalnost minimalnog pojacanja

Fastep = 0.01*Fae;
Fbstep = 0.01*Fbe;

tu1_step1 = 1000;
%tu1_step2 = 6000;
tu1_step2 = 10000; %za min pojacanje
%tu1_step2 = 8000; %za poremecaj

tu2_step1 = 1000;
%tu2_step2 = 6000;
tu2_step2 = 10000; %za min pojacanje
%tu2_step2 = 8000; %za poremecaj
%% Skaliranje

pHmax =  0.0596; 
%pHmax =  0.001*0.0596; %za direkcionalnost minimalnog pojacanja
xamax = 2.0819e-08;
%xamax = pHstep*0.0003816; %za direkcionalnost minimalnog pojacanja
Famax = 0.1*Fae;
Fbmax = 0.1*Fbe;

Dy = [pHmax 0; 0 xamax];
Du = [Famax 0; 0 Fbmax];

%% Decentralizovan kontroler: zanemarivanje sprežnih elemenata
s = tf('s');

skaliran_Kstar_dec_z = zpk([0.80709*(s+0.002)*(s+0.001667)/s/(s+0.00181) 0; 0 0.95559*(s+0.001667)/s]);

C = [0 1; 1 0];

Kstar_dec_z = Du*skaliran_Kstar_dec_z*inv(Dy);
K_dec_z   = inv(C)*Kstar_dec_z;
%% Decentralizovan kontroler: statičko dekuplovanje

skaliran_K_dec_stat = [3.9687*(s+0.002)/s 3.9687*(s+0.002)/s; 
    4.9117*(s+0.001667)/s 3.9423*(s+0.001667)/s];
K_dec_stat = Du*skaliran_K_dec_stat*inv(Dy);
%% Decentralizovan kontroler: dekuplovanje na učestanosti propusnog opsega

skaliran_K_dec_w0 = [3.9687*(s+0.002)/s 3.9687*(s+0.002)/s;
    4.7139*(s+0.001667)/s 3.7445*(s+0.001667)/s];
K_dec_w0 = Du*skaliran_K_dec_w0*inv(Dy);

%% Hinf kontroler

skaliran_K_hinf = [(5.269*s^2 + 0.5178*s + 0.001015)/(s^3 + 1.216*s^2 + 0.1131*s + 1.923e-09) (4.697*s^2 + 0.5161*s + 0.001014)/(s^3 + 1.216*s^2 + 0.1131*s + 1.923e-09);
    (5.524*s^2 + 0.5966*s + 0.00106)/(s^3 + 1.216*s^2 + 0.1131*s + 1.923e-09) (4.818*s^2 + 0.4692*s + 0.0008497)/(s^3 + 1.216*s^2 + 0.1131*s + 1.923e-09)];
K_hinf = Du*skaliran_K_hinf*inv(Dy);

%% Kontroler sa potpunom inverzijom dinamike
skaliran_K_inv = [3.7664*(s+0.002)/s 3.7664*(s+0.002)/s; 4.4594*(s+0.00181)/s 3.5038*(s+0.001848)/s];
K_inv = Du*skaliran_K_inv*inv(Dy);
%% Pokretanje simulacije
sim('kontrol_nelin_svi.slx');

%% nelinearni sistem

% min_xa_out_plot = 0;% min(xa_out);
% max_xa_out_plot = max(xa_out);
% delta_xa_out = max_xa_out_plot-xa_out(1);
% 
% min_pH_out_plot = 0; %min(pH_out);
% max_pH_out_plot = max(pH_out);
% delta_pH_out = max_pH_out_plot-pH_out(1);
% 
% min_Fa_out_plot = 0; %min(Fa_out);
% max_Fa_out_plot = max(Fa_out);
% 
% min_Fb_out_plot = 0; % min(Fb_out);
% max_Fb_out_plot = max(Fb_out);
%%
figure;
hold all;
xlim([0 2e4])
plot(t_out,pH_ref_out,'k--','LineWidth',1)
plot(t_out,pH_dec_z_out,'r','LineWidth',1)
plot(t_out,pH_dec_stat_out,'g','LineWidth',1)
plot(t_out,pH_dec_w0_out,'b','LineWidth',1)
plot(t_out,pH_hinf_out,'c','LineWidth',1)
plot(t_out,pH_inv_out,'m','LineWidth',1)

grid
xlabel('vreme [s]')
ylabel('pH [au]')
title('pH vrednost tečnosti u reaktoru')
legend('referentna vrednost','dec, zan','dec stat','dec w0','hinf','inv')
%%
figure;
hold all;
xlim([0 2e4])
plot(t_out,xa_ref_out,'k--','LineWidth',1)
plot(t_out,xa_dec_z_out,'r','LineWidth',1)
plot(t_out,xa_dec_stat_out,'g','LineWidth',1)
plot(t_out,xa_dec_w0_out,'b','LineWidth',1)
plot(t_out,xa_hinf_out,'c','LineWidth',1)
plot(t_out,xa_inv_out,'m','LineWidth',1)

grid
xlabel('vreme [s]')
ylabel('xa [mol/l]')
title('Koncentracija kiseline u reaktoru')
legend('referentna vrednost','dec, zan','dec stat','dec w0','hinf','inv')
%%
figure;
hold all;
xlim([0 2e4])
plot(t_out,Fa_dec_z_out,'r','LineWidth',1)
plot(t_out,Fa_dec_stat_out,'g','LineWidth',1)
plot(t_out,Fa_dec_w0_out,'b','LineWidth',1)
plot(t_out,Fa_hinf_out,'c','LineWidth',1)
plot(t_out,Fa_inv_out,'m','LineWidth',1)

grid
ylabel('Fa [l/s]')
xlabel('vreme [s]')
title('Protok kiseline u reaktoru')
legend('dec, zan','dec stat','dec w0','hinf','inv')
%%
figure;
hold all;
xlim([0 2e4])
plot(t_out,Fb_dec_z_out,'r','LineWidth',1)
plot(t_out,Fb_dec_stat_out,'g','LineWidth',1)
plot(t_out,Fb_dec_w0_out,'b','LineWidth',1)
plot(t_out,Fb_hinf_out,'c','LineWidth',1)
plot(t_out,Fb_inv_out,'m','LineWidth',1)

grid
ylabel('Fb [l/s]')
xlabel('vreme [s]')
title('Protok baze u reaktoru')
legend('dec, zan','dec stat','dec w0','hinf','inv')