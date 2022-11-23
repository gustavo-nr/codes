%% Seguidor de referência variante e rejeitador de distúrbio

clc
clear
close all

addpath('Matrizes\')
addpath("Imagens\Controle Moderno\")

%% Importa matrizes
A = importdata('matrix_A1lin.txt');
B2 = importdata('matrix_B1lin.txt');
B1 = -B2;
C = importdata('matrix_C1.txt');
D = importdata('matrix_D1.txt');

%% Define a referência constante
x_r = [pi/3 pi/4 pi/5 0 0 0];
seq = 'ZXZ';
q_r = angle2quat(x_r(1),x_r(2),x_r(3),seq);
xq_r = [q_r(2:4),x_r(4:6)];

%% Regulador por alocação de polos
p1 = -0.3+0.3i;
p2 = -0.4+0.3i;
p = [p1 conj(p1) p2 conj(p2) -0.5 -0.4];
K = place(A,B2,p);
F = A-B2*K;
T=40;

%% Matrizes de referência e disturbios
A_r = [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 0 0];

A_w = [-41.9938648159790	-35.9947506215021	-5795.14933552378;
233304.083573796	199975.061457068	32195941.2600482;
-1448.78938342553	-1241.82029553629	-199932.796599317];

%% Desenvolvimento da malha fechada
A_z = zeros(length(A),length(A));

A_o = [A_w zeros(length(A_w),length(A_r));...
    zeros(length(A_r),length(A_w)) A_r];

F2 = [B1 (A-A_r)];
C_bar = [eye(3),zeros(3,3)];
K_ex = inv(C_bar*inv(F)*B2)*C_bar*inv(F)*F2;

A_y = [B1 B2*K]-B2*K_ex;

%% Definição do sistema
A_bar = [F A_y;zeros(9,6) A_o];
x_bar_0 = [0 0 0 0 0 0 ...
    0 0 0 ...
    xq_r(1) xq_r(2) xq_r(3) 0 0 0]';
sys = ss(A_bar,x_bar_0,A_bar,x_bar_0);

%% Simulação do sistema
[y,t]=step(sys,60);

q0 = q_r(1);
q1 = y(:,1);
q2 = y(:,2);
q3 = y(:,3);
w1 = y(:,4);
w2 = y(:,5);
w3 = y(:,6);

[psi,theta,phi] = quat2eulang(q0,q1,q2,q3,seq);

figure(1)
plot(t,w1,LineWidth=1.20)
hold on
plot(t,w2,LineWidth=1.20)
plot(t,w3,LineWidth=1.20)
grid on
title("Velocidade angular ao longo do tempo")
xlabel("Tempo [s]")
ylabel("Velocidade angular [rad/s]")
legend('w1','w2','w3')
hold off

baseFileName = sprintf('Image_%s.png', "seg_var_w_t");
fullFileName = fullfile("Imagens\Controle Moderno\", baseFileName);
saveas(1, fullFileName);

figure(2)
plot(t,q1,LineWidth=1.20)
hold on
plot(t,q2,LineWidth=1.20)
plot(t,q3,LineWidth=1.20)
grid on
title("Posição angular em quatérnios ao longo do tempo")
xlabel("Tempo [s]")
ylabel("Posição angular [rad]")
legend('q1','q2','q3','Location','northwest')
hold off

baseFileName = sprintf('Image_%s.png', "seg_var_quat_t");
fullFileName = fullfile("Imagens\Controle Moderno\", baseFileName);
saveas(2, fullFileName);

figure(3)
plot(t,psi,LineWidth=1.20)
hold on
plot(t,theta,LineWidth=1.20)
plot(t,phi,LineWidth=1.20)
grid on
title("Posição angular em ângulos de Euler ao longo do tempoo")
xlabel("Tempo [s]")
ylabel("Posição angular [rad]")
legend('psi','theta','phi','Location','southeast')
hold off

baseFileName = sprintf('Image_%s.png', "seg_var_euler_t");
fullFileName = fullfile("Imagens\Controle Moderno", baseFileName);
saveas(3, fullFileName);



