clear all
close all 
clc
 
%% Constantes
N = 10;
fe = 1;
Te = 1/fe;
time = (Te:Te:N);

delta_T = 1;
vit_remplissage = 0.1;

%% Modele dynamique
L = 1;
%Z_th = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]';
%Z = [0.0, 0.09, 0.21, 0.35, 0.39, 0.54, 0.62, 0.74, 0.82, 0.88]';

Z_th = vit_remplissage * time';
std_V_reel = 0.1;
V_reel = std_V_reel * randn(N,1);
Z = Z_th + V_reel;

dimEtat = 2;
A = [ 1 delta_T; 
      0 1];
B = 0;
u = 0;

%Covariance du bruit d'etat
Q = 10*eye(dimEtat);

%Covariance du bruit de mesure
dimMesure = 1;
H = [1 0];
R = 1E-3*eye(dimMesure);

%Estimation initiale
X0_estime = [0. .1]';
P0_estime = [1E-3 0; 0 1E-3];

%Creation des vecteurs
X_kalman = zeros(N+1,dimEtat);

Innovation = zeros(N+1,1);
K = zeros(N+1,2);

P = zeros(N+1,2,2);

%Initialisation
P_courant = P0_estime;
X_courant = X0_estime;

P(1,:,:) = P0_estime;
X_kalman(1,:) = X0_estime;
K(1,:) = (P_courant*H') / ((H*P_courant)*H' + R);

for k = 1:N*fe
    %Prediction
    X_courant = A*X_courant + B*u;
    P_courant = A*P_courant*A' + Q;
    
    %Correction
    K_courant = (P_courant*H')/(H*P_courant*H' + R);
    Innovation_courant = Z(k) - H*X_courant;
    X_courant = X_courant + K_courant*Innovation_courant;
    P_courant = (eye(dimEtat) - K_courant*H)*P_courant;
    
    %Sauvegarde des donnees
    X_kalman(k+1,:) = X_courant';
    K(k+1,:) = K_courant';
    Innovation(k+1,:) = Innovation_courant';
    P(k+1,:,:) = P_courant;
end

figure
subplot(411)
plot(time,K(2:end,2))
title('Gain du filtre de Kalman')
grid on

subplot(412)
plot(time,P(2:end,1))
title('Variance de l''erreur d''estimation du niveau d''eau (P11)')
grid on

subplot(413)
plot(time,Innovation(2:end))
title('Innovation')
grid on

subplot(414)
plot(time,P(2:end,2))
title('Variance de l''erreur d''estimation de la vitesse (P22)')
grid on


figure
subplot(211)
plot(time,Z_th,'r')
hold on
plot(time,Z,'ob')
hold on
plot(time, X_kalman(2:end,1),'g')
title('Niveau d''eau')
grid on
legend('Mesures theoriques', 'Mesures reelles', 'Estimees')

subplot(212)
plot(time,vit_remplissage*ones(N,1),'r')
hold on
plot(time,X_kalman(2:end,2),'g')
title('Vitesse de remplissage')
grid on
legend('Mesures theoriques', 'Estimees')
