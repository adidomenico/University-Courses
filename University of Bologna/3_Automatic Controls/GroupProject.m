clear all
close all
clc
% Progetto di Controlli Automatici

%Yuri Noviello 0000873239   
%Roberto Calvani 0000881408
%Alessio Di Domenico 0000882684
%Tommaso Ferrari 0000873705

% Traccia 3C, Gruppo Y

%% Parametri del sistema
mi = 0.15; %massa[kg]
ei = 0.47; %distanza dalla massa[m] (+-0.1)
ei_m = 0.37; %ei-0.1
ei_p= 0.57; %ei+0.1
beta = 0.3; %coefficiente moltiplicativo
Ie = 0.75; %momento di inerzia[kg*m*m]
mud = 0.012; %coefficiente attrito dinamico

we = 1590; %velocita di equilibrio di equilibrio [rad/s]

%Ricavo u_e (coppia applicata al rotore di equilibrio)
tau_e = mud*mi*ei^2*we^2 + beta*we; 

%equazioni che descrivono il sistema
%f1=w
%f2=-beta/(miei^2+Ie)*w - (mud*mi*e^2)/(miei^2+Ie)*w^2 + tau/(miei^2+Ie)

% Matrici di linearizzazione del sistema dinamico
AA = [0, 1; 0, -beta/(mi*ei^2+Ie) - 2*(mud*mi*ei^2)/(mi*ei^2+Ie)*we ];
BB = [0; 1/(mi*ei^2+Ie)];
CC = [0, 1];
DD = 0; 

%% Funzione di Trasferimento equivalente
s = tf('s');
GG = (1/(mi*ei^2+Ie)) / (s + beta/(mi*ei^2+Ie) + 2*(mud*mi*ei^2)/(mi*ei^2+Ie)*we) ;
GG_m = (1/(mi*ei_m^2+Ie)) / (s + beta/(mi*ei_m^2+Ie) + 2*(mud*mi*ei_m^2)/(mi*ei_m^2+Ie)*we) ;
GG_p =  (1/(mi*ei_p^2+Ie)) / (s + beta/(mi*ei_p^2+Ie) + 2*(mud*mi*ei_p^2)/(mi*ei_p^2+Ie)*we) ;

figure(1);
pzmap(GG);
grid on, zoom on
title ("ZERI E POLI") %figura che rappresenta la posizione degli zeri e poli

figure(2);
bode(GG);
grid on, zoom on
title ("BODE G(s)") %figura che rappresenta diagramma di bode della funzione

%% Specifiche di progetto
W=4; % [m] costante moltiplicativa del gradino
S_100= 0.01; % Sovraelongazione percentuale

%ricavo csi da S_100
xi= sqrt(log(S_100)^2/(pi^2+log(S_100)^2)); %csi: coefficiente di smorzamento

Mf=xi*100; %Margine di fase
display(Mf)

T_a1Max=4; %[s] Tempo di assestamento all' 1_perc

%Ricavo omega_cMin per il margine di Ampiezza
omega_cMin=4.6/(T_a1Max*xi);
display(omega_cMin);

%Specifiche sul rumore di misura
omega_n=120;
A_n=-20*log10(30); %-30

%% Rappresentazione grafica della risposta al gradino in relazione ai vincoli di progetto
figure(3); 
T_simulation=20; %tempo di simulazione
step(GG,T_simulation,'b')
hold on
title ("Risposta al Gradino") 
% Sovraelongazione
patch([0,T_simulation,T_simulation,0],[W*(1+S_100),W*(1+S_100),W+1,W+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on; ylim([0,W+1]);
% Tempo di Assestamento
patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1-0.01),W*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1+0.01),W*(1+0.01),W+1,W+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

hold off
Legend=["G(s) Risposta al gradino";"Vincolo Sovraelongazione"; "Vincolo Tempo di Assestamento"];
legend(Legend);


%% Rappresentazione grafica dei digarammi di Bode in relazione ai vincoli di progetto
omega_plot_min=10^(-1);
omega_plot_max=10^4;

omega_cMax=omega_n;

[Mag,phase,omega]=bode(GG,{omega_plot_min,omega_plot_max},'k');
figure(4)
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Diagramma di Bode G(jw)")

%% Realizzo la Funzione di Trasferimento estesa , aggiungendo un polo per avere errore a regime nullo (Regolazione statica)
GG_e=GG/s; % Ge= G(s) * Rs(s), con Rs(s) = 1/s 
GG_e_p=GG_p/s;
GG_e_m=GG_m/s;
[Mag,phase,omega]=bode(GG_e,{omega_plot_min,omega_plot_max},'k');
figure(5)
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G_e(j\omega)" ];

legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Diagramma di Bode Ge(jw)")

%% Realizzazione della rete anticipatrice
omega_cStar= (omega_cMin+44); %Scelgo omega_cStar
argGe=-(atan(2*omega_cStar / -omega_cStar^2)+180);
phi_star= (Mf+3) - 180 -(argGe);display(phi_star); %Calcolo phi_star
display(omega_cStar);

%M_star=10 ^(-|Ge(jw_c*)|db / 20 ) ;
M_star=30;

display(1/M_star)
display (cosd(phi_star))
%Mi assicuro che la disequazione cos(phi_star)>1/M_star sia rispettata

tau= (M_star - cosd(phi_star))/(omega_cStar * sind(phi_star)); 
atau=(cosd(phi_star) - 1/M_star)/(omega_cStar * sind(phi_star));

R_ant=(1+tau*s)/(1+atau*s);
display(R_ant);
figure(6)
bode(R_ant)
title("Rete anticipatrice")

%% Aggiungo la rete anticipatrice alla G_estesa

[Mag,phase,omega]=bode(GG_e*R_ant,{omega_plot_min,omega_plot_max},'k');
figure(7)
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G_e R_{ant}(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e R_{ant}(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Diagramma di Bode Ge(jw)R_ant(jw)")
zoom on

% Controlliamo la risposta al gradino
figure(8)
L=GG/s*R_ant;
display(L);
step(L/(1+L)) %Funzione di sensivita complementare

%% Aggiungiamo il gain e il polo ad alta frequenza
%Control System Design,
gain = 4;
R2 = 115/(s+115);

L     =  GG_e*R_ant*gain*R2;
L_p =  GG_e_p*R_ant*gain*R2; 
L_m = GG_e_m*R_ant*gain*R2; 

%Grafichiamo L(jw), GG_e(jw), GG(jw)
[Mag,phase,omega]=bode(L,{omega_plot_min,omega_plot_max},'k');
[Mag1,phase1,omega1]=bode(GG_e,{omega_plot_min,omega_plot_max},'r');
[Mag2,phase2,omega2]=bode(GG,{omega_plot_min,omega_plot_max},'k');

figure(9)
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "L(j\omega)"
    "GG_e(j\omega)"
    "GG(j\omega)"];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);
margin(Mag1,phase1,omega1);
margin(Mag2,phase2,omega2);
grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["L(j\omega)" ;"GG_e(j\omega)" ;"GG(j\omega)" ;"M_f"] ;
legend(Legend_arg)
hold off
title("Diagramma di Bode di L(jw), G_e(jw), G(jw) ")
zoom on

figure(10)
step(L/(1+L),10) %Funzione di sensivita complementare
grid on
zoom on 

%% Risposta al Gradino

figure(11)
[y_step,t_step]=step(W*L/(1+L),20);
[y_step1,t_step1]=step(W*L_p/(1+L_p),20);
[y_step2,t_step2]=step(W*L_m/(1+L_m),20);
plot(t_step,y_step)
plot(t_step1,y_step1)
plot(t_step2,y_step2)
plot(t_step,y_step,'r', t_step1,y_step1,'k', t_step2,y_step2,'b')

patch([0,4,4,0],[W*(1+S_100),W*(1+S_100),W+1,W+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on; ylim([0,W+1]);

patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1-0.01),W*(1-0.01),0,0],'p','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1+0.01),W*(1+0.01),W+1,W+1],'p','FaceAlpha',0.1,'EdgeAlpha',0.1);

hold off
Legend=["F(s) Risposta al Gradino" ;"F(s) +0.1 Risposta al Gradino";"F(s) -0.1 Risposta al Gradino" ;"Vincolo di sovraelongazione"; "Vincolo del Tempo di Assestamento"];
legend(Legend,'location','SouthEast');

zoom on


%% Definisco i dati per le simulazione con Simulink

R=(1/s)*R_ant*gain*R2; 

[n_r,d_r]=tfdata(R);
num_r=n_r{1};
den_r=d_r{1};

W=4;

x0 = [-90, 0];
xe = [-90, we];
ue=tau_e;
ye= we;

open('rotore.slx');
sim('rotore.slx', 20);
%% Punto opzionale 1
%Valori soglia minima e massima
%xe = [-90, we-10];
%xe= [-90 , we+5];

%% Punto opzionale 2
%Soglia massima per W
%W=6;