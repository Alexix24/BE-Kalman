clear all
close all 
clc
 
%% Constantes
N = 10;
fe = 1;
Te = 1/fe;
time = (Te:Te:N);

%% Modele statique
L = 1;
Z_th = L * ones(N,1);

%Z = [0.9, 0.8, 1.1, 1, 0.95, 1.05, 1.2, 0.9, 0.85, 1.15]';

std_V_reel = 0.1;
V_reel = std_V_reel * randn(N,1);
Z = Z_th + V_reel;

dimEtat = 1;
A = 1;  %Matrice de transition
B = 0;

u = 0;  %Pas de commande dans ce système

%Covariance du bruit de mesure
dimMesure = 1;
H = 1;
R = 100*eye(dimMesure);

%Covariance du bruit d'etat
Q = 1E-3*eye(dimEtat);



%Estimation initiale
X0_estime = 0.5;
P0_estime = 1E-3;

%Creation des vecteurs
X_kalman = zeros(N+1,dimEtat);

Innovation = zeros(N+1,1);
K = zeros(N+1,1);

P = zeros(N+1,1);
P(1) = P0_estime;

X_kalman(1) = X0_estime;

K(1) = (P(1)*H') / (H*P(1)*H' + R);

%Algorithme de Kalman
for k = 1:N*fe
    
    X_kalman(k+1) = A*X_kalman(k) + B*u;
    P(k+1) = A*P(k)*A' + Q;
    
    K(k+1) = (P(k+1)*H')/(H*P(k+1)*H' + R);
    Innovation(k+1) = Z(k) - H*X_kalman(k+1);
    X_kalman(k+1) = X_kalman(k+1) + K(k+1)*(Innovation(k+1));
    P(k+1) = (eye(dimEtat) - K(k+1)*H)*P(k+1);
    
end

figure
subplot(311)
plot(time,K(2:end))
title('Gain du filtre de Kalman')
grid on

subplot(312)
plot(time,P(2:end))
title('Variance de l''erreur d estimation(P)')
grid on

subplot(313)
plot(time,Innovation(2:end))
title('Innovation')
grid on

figure
plot(time,Z_th,'r')
hold on
plot(time,Z,'ob')
hold on
plot(time, X_kalman(2:end),'g')
title('Niveau d''eau')
grid on
legend('mesure theorique', 'mesure reelle', 'estimees')